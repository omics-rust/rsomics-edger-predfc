//! edgeR predFC: predictive log2-fold-change coefficients from a prior-count
//! augmented negative-binomial GLM fit. Method: McCarthy, Chen & Smyth (2012),
//! NAR 40:4288-4297. A library-size-scaled prior count is added to every count
//! and folded into the offset, the NB-GLM (log link) is fit to convergence by
//! Levenberg-damped Fisher scoring, and its coefficients are returned in log2.

use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

const LN2: f64 = std::f64::consts::LN_2;
const MAXIT: usize = 200;
const TOL: f64 = 1e-10;

pub struct Matrix {
    pub gene_col: String,
    pub sample_names: Vec<String>,
    pub genes: Vec<String>,
    pub counts: Vec<f64>,
    pub n_samples: usize,
}

impl Matrix {
    pub fn load(path: &Path) -> Result<Self> {
        let file = File::open(path)
            .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?;
        let mut lines = BufReader::new(file).lines();
        let header = lines
            .next()
            .ok_or_else(|| RsomicsError::InvalidInput("empty count matrix".into()))?
            .map_err(RsomicsError::Io)?;
        let mut hcols = header.split('\t');
        let gene_col = hcols.next().unwrap_or("gene").to_string();
        let sample_names: Vec<String> = hcols.map(str::to_string).collect();
        let n_samples = sample_names.len();
        if n_samples == 0 {
            return Err(RsomicsError::InvalidInput(
                "count matrix has no sample columns".into(),
            ));
        }
        let mut genes = Vec::new();
        let mut counts = Vec::new();
        for line in lines {
            let line = line.map_err(RsomicsError::Io)?;
            if line.is_empty() {
                continue;
            }
            let mut fields = line.split('\t');
            let gene = fields
                .next()
                .ok_or_else(|| RsomicsError::InvalidInput("row without a gene id".into()))?;
            genes.push(gene.to_string());
            let before = counts.len();
            for f in fields {
                counts.push(f.parse::<f64>().map_err(|_| {
                    RsomicsError::InvalidInput(format!("non-numeric count '{f}' for gene {gene}"))
                })?);
            }
            if counts.len() - before != n_samples {
                return Err(RsomicsError::InvalidInput(format!(
                    "gene {gene}: {} values, header has {n_samples} samples",
                    counts.len() - before
                )));
            }
        }
        Ok(Self {
            gene_col,
            sample_names,
            genes,
            counts,
            n_samples,
        })
    }

    pub fn n_genes(&self) -> usize {
        self.genes.len()
    }
    fn row(&self, g: usize) -> &[f64] {
        &self.counts[g * self.n_samples..(g + 1) * self.n_samples]
    }
}

/// Design matrix: row-major n_samples × n_coef, plus column names from the header.
pub struct Design {
    pub data: Vec<f64>,
    pub n_samples: usize,
    pub n_coef: usize,
    pub coef_names: Vec<String>,
}

impl Design {
    fn load(path: &Path) -> Result<Self> {
        let file = File::open(path)
            .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?;
        let mut lines = BufReader::new(file).lines();
        let header = lines
            .next()
            .ok_or_else(|| RsomicsError::InvalidInput("empty design matrix".into()))?
            .map_err(RsomicsError::Io)?;
        let coef_names: Vec<String> = header.split('\t').map(str::to_string).collect();
        let n_coef = coef_names.len();
        let mut data = Vec::new();
        let mut n_samples = 0;
        for line in lines {
            let line = line.map_err(RsomicsError::Io)?;
            if line.is_empty() {
                continue;
            }
            let before = data.len();
            for f in line.split('\t') {
                data.push(f.parse::<f64>().map_err(|_| {
                    RsomicsError::InvalidInput(format!("non-numeric design value '{f}'"))
                })?);
            }
            if data.len() - before != n_coef {
                return Err(RsomicsError::InvalidInput(format!(
                    "design row {n_samples}: {} values, header has {n_coef} columns",
                    data.len() - before
                )));
            }
            n_samples += 1;
        }
        Ok(Self {
            data,
            n_samples,
            n_coef,
            coef_names,
        })
    }

    fn row(&self, s: usize) -> &[f64] {
        &self.data[s * self.n_coef..(s + 1) * self.n_coef]
    }
}

/// One value per row, taking the last tab-separated column so a `gene<TAB>value`
/// or a bare-value file both parse. Comment and blank lines are skipped.
fn load_column(path: &Path, expected: usize, what: &str) -> Result<Vec<f64>> {
    let file = File::open(path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?;
    let mut v = Vec::with_capacity(expected);
    for line in BufReader::new(file).lines() {
        let line = line.map_err(RsomicsError::Io)?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let val = line.rsplit('\t').next().unwrap_or(line);
        v.push(
            val.parse::<f64>()
                .map_err(|_| RsomicsError::InvalidInput(format!("non-numeric {what} '{val}'")))?,
        );
    }
    if v.len() != expected {
        return Err(RsomicsError::InvalidInput(format!(
            "{} {what} values, expected {expected}",
            v.len()
        )));
    }
    Ok(v)
}

/// Reusable per-thread scratch so the hot IRLS loop allocates nothing per gene.
struct Scratch {
    beta: Vec<f64>,
    trial: Vec<f64>,
    mu: Vec<f64>,
    mu_t: Vec<f64>,
    xtwx: Vec<f64>,
    a: Vec<f64>,
    xtr: Vec<f64>,
    rhs: Vec<f64>,
    step: Vec<f64>,
    row_aug: Vec<f64>,
}

impl Scratch {
    fn new(n: usize, p: usize) -> Self {
        Scratch {
            beta: vec![0.0; p],
            trial: vec![0.0; p],
            mu: vec![0.0; n],
            mu_t: vec![0.0; n],
            xtwx: vec![0.0; p * p],
            a: vec![0.0; p * p],
            xtr: vec![0.0; p],
            rhs: vec![0.0; p],
            step: vec![0.0; p],
            row_aug: vec![0.0; n],
        }
    }
}

/// NB deviance (edgeR nbinomDeviance), the IRLS line-search objective.
fn nb_deviance(y: &[f64], mu: &[f64], dispersion: f64) -> f64 {
    let r = if dispersion > 0.0 {
        1.0 / dispersion
    } else {
        f64::INFINITY
    };
    let mut dev = 0.0;
    for (&yi, &mui) in y.iter().zip(mu) {
        let term_y = if yi > 0.0 { yi * (yi / mui).ln() } else { 0.0 };
        dev += if dispersion > 0.0 {
            term_y - (yi + r) * ((yi + r) / (mui + r)).ln()
        } else {
            term_y - (yi - mui)
        };
    }
    2.0 * dev
}

fn eta_mu(design: &Design, offset: &[f64], beta: &[f64], mu: &mut [f64]) {
    for s in 0..design.n_samples {
        let mut eta = offset[s];
        for (&x, &b) in design.row(s).iter().zip(beta) {
            eta += x * b;
        }
        mu[s] = eta.exp();
    }
}

/// One-group NB fit (edgeR mglmOneGroup): the single coefficient β with
/// μ[j] = exp(β + offset[j]), on the natural-log scale. Seeds the per-gene start.
fn mglm_one_group(row: &[f64], offset: &[f64], dispersion: f64) -> f64 {
    let total: f64 = row.iter().sum();
    if total == 0.0 {
        return f64::NEG_INFINITY;
    }
    let mean_off = offset.iter().sum::<f64>() / offset.len() as f64;
    let mut beta = (total / row.len() as f64).ln() - mean_off;
    for _ in 0..MAXIT {
        let mut dl = 0.0;
        let mut info = 0.0;
        for (&y, &off) in row.iter().zip(offset) {
            let mu = (beta + off).exp();
            let denom = 1.0 + mu * dispersion;
            dl += (y - mu) / denom;
            info += mu / denom;
        }
        let s = dl / info;
        beta += s;
        if s.abs() < TOL {
            break;
        }
    }
    beta
}

/// Per-design start direction d = (XᵀX)⁻¹Xᵀ1, so edgeR's null-method start for a
/// gene is b0·d (the OLS projection of a constant linear predictor onto X).
fn start_direction(design: &Design) -> Vec<f64> {
    let p = design.n_coef;
    let mut xtx = vec![0.0f64; p * p];
    let mut xt1 = vec![0.0f64; p];
    for s in 0..design.n_samples {
        let xr = design.row(s);
        for (j, &xj) in xr.iter().enumerate() {
            xt1[j] += xj;
            let rowj = &mut xtx[j * p..j * p + p];
            for (rk, &xk) in rowj.iter_mut().zip(xr) {
                *rk += xj * xk;
            }
        }
    }
    let mut d = vec![0.0f64; p];
    if !solve(&mut xtx, &mut xt1, &mut d, p) {
        return vec![0.0f64; p];
    }
    d
}

/// IRLS NB-GLM (log link, per-sample offset), Fisher scoring with a Levenberg
/// ridge — edgeR mglmLevenberg, run tight to the β-MLE. Returns the coefficients.
fn fit_nb_glm(
    y: &[f64],
    design: &Design,
    offset: &[f64],
    dispersion: f64,
    start_dir: &[f64],
    sc: &mut Scratch,
) -> Vec<f64> {
    let n = design.n_samples;
    let p = design.n_coef;

    let b0 = mglm_one_group(y, offset, dispersion);
    for (b, &d) in sc.beta[..p].iter_mut().zip(start_dir) {
        *b = b0 * d;
    }

    eta_mu(design, offset, &sc.beta[..p], &mut sc.mu);
    let mut dev = nb_deviance(y, &sc.mu[..n], dispersion);
    let mut lambda = 0.0f64;

    for _ in 0..MAXIT {
        for v in sc.xtwx[..p * p].iter_mut() {
            *v = 0.0;
        }
        for v in sc.xtr[..p].iter_mut() {
            *v = 0.0;
        }
        let rows = design.data[..n * p].chunks_exact(p);
        for ((&yi, &mui), xr) in y.iter().zip(&sc.mu[..n]).zip(rows) {
            let denom = 1.0 + dispersion * mui;
            let w = mui / denom;
            let resid = (yi - mui) / denom;
            for (j, &xj) in xr.iter().enumerate() {
                sc.xtr[j] += xj * resid;
                let xjw = xj * w;
                let rowj = &mut sc.xtwx[j * p..j * p + p];
                for (rk, &xk) in rowj.iter_mut().zip(xr) {
                    *rk += xjw * xk;
                }
            }
        }
        let mut accepted = false;
        for _ in 0..20 {
            sc.a[..p * p].copy_from_slice(&sc.xtwx[..p * p]);
            for d in 0..p {
                sc.a[d * p + d] += lambda * sc.xtwx[d * p + d].max(1e-6);
            }
            sc.rhs[..p].copy_from_slice(&sc.xtr[..p]);
            if !solve(&mut sc.a[..p * p], &mut sc.rhs[..p], &mut sc.step[..p], p) {
                lambda = if lambda == 0.0 { 1.0 } else { lambda * 2.0 };
                continue;
            }
            for (t, (&b, &s)) in sc.trial[..p]
                .iter_mut()
                .zip(sc.beta[..p].iter().zip(&sc.step[..p]))
            {
                *t = b + s;
            }
            eta_mu(design, offset, &sc.trial[..p], &mut sc.mu_t);
            let dev_t = nb_deviance(y, &sc.mu_t[..n], dispersion);
            if dev_t <= dev + 1e-8 * (1.0 + dev.abs()) {
                let max_step = sc.step[..p].iter().fold(0.0f64, |m, s| m.max(s.abs()));
                sc.beta[..p].copy_from_slice(&sc.trial[..p]);
                sc.mu[..n].copy_from_slice(&sc.mu_t[..n]);
                dev = dev_t;
                lambda *= 0.5;
                accepted = true;
                if max_step < TOL {
                    return sc.beta[..p].to_vec();
                }
                break;
            }
            lambda = if lambda == 0.0 { 1.0 } else { lambda * 4.0 };
        }
        if !accepted {
            break;
        }
    }
    sc.beta[..p].to_vec()
}

/// Solve A x = b (A row-major p×p) by Gaussian elimination with partial pivoting.
/// `a` and `rhs` are clobbered. Returns false if singular.
fn solve(a: &mut [f64], rhs: &mut [f64], x: &mut [f64], p: usize) -> bool {
    for col in 0..p {
        let mut piv = col;
        let mut best = a[col * p + col].abs();
        for r in (col + 1)..p {
            let v = a[r * p + col].abs();
            if v > best {
                best = v;
                piv = r;
            }
        }
        if best < 1e-12 {
            return false;
        }
        if piv != col {
            for k in 0..p {
                a.swap(col * p + k, piv * p + k);
            }
            rhs.swap(col, piv);
        }
        let d = a[col * p + col];
        for r in (col + 1)..p {
            let f = a[r * p + col] / d;
            if f == 0.0 {
                continue;
            }
            for k in col..p {
                a[r * p + k] -= f * a[col * p + k];
            }
            rhs[r] -= f * rhs[col];
        }
    }
    for col in (0..p).rev() {
        let mut s = rhs[col];
        for k in (col + 1)..p {
            s -= a[col * p + k] * x[k];
        }
        x[col] = s / a[col * p + col];
    }
    true
}

pub struct PredFcArgs<'a> {
    pub counts: &'a Path,
    pub design: &'a Path,
    pub dispersion: f64,
    pub dispersion_file: Option<&'a Path>,
    pub norm_factors: Option<&'a Path>,
    pub offset_file: Option<&'a Path>,
    pub prior_count: f64,
}

pub fn predfc(args: &PredFcArgs, output: &mut dyn Write) -> Result<u64> {
    let m = Matrix::load(args.counts)?;
    let design = Design::load(args.design)?;
    if design.n_samples != m.n_samples {
        return Err(RsomicsError::InvalidInput(format!(
            "design has {} rows but matrix has {} samples",
            design.n_samples, m.n_samples
        )));
    }

    let mut lib = vec![0.0f64; m.n_samples];
    for row in m.counts.chunks_exact(m.n_samples) {
        for (s, &c) in lib.iter_mut().zip(row) {
            *s += c;
        }
    }

    // The fit runs against an offset of log(effective library size). edgeR lets
    // it come straight from a DGEList (log(lib·normfactors)) or be supplied; an
    // explicit offset file wins, else lib·normfactors, else the raw lib sizes.
    let offset: Vec<f64> = if let Some(p) = args.offset_file {
        load_column(p, m.n_samples, "offset")?
    } else {
        let nf = match args.norm_factors {
            Some(p) => load_column(p, m.n_samples, "norm factor")?,
            None => vec![1.0; m.n_samples],
        };
        lib.iter().zip(&nf).map(|(&l, &f)| (l * f).ln()).collect()
    };

    // A zero effective library size makes offset = log(0) = -Inf, which poisons
    // the GLM into an all-NaN fit. edgeR bails here (.compressOffsets: "offsets
    // must be finite values"); so do we, rather than ship a NaN table with a
    // success exit. Covers the --offset override path too (a supplied ±Inf/NaN).
    if let Some(s) = offset.iter().position(|o| !o.is_finite()) {
        let id = m
            .sample_names
            .get(s)
            .map_or_else(|| s.to_string(), Clone::clone);
        return Err(RsomicsError::InvalidInput(format!(
            "zero library size for sample {id}: offsets must be finite values"
        )));
    }

    let eff_lib: Vec<f64> = offset.iter().map(|&o| o.exp()).collect();

    let mean_eff = eff_lib.iter().sum::<f64>() / eff_lib.len() as f64;
    let prior: Vec<f64> = eff_lib
        .iter()
        .map(|&l| args.prior_count * l / mean_eff)
        .collect();
    let offset_aug: Vec<f64> = eff_lib
        .iter()
        .zip(&prior)
        .map(|(&l, &p)| (l + 2.0 * p).ln())
        .collect();

    let dispersions = match args.dispersion_file {
        Some(p) => load_column(p, m.n_genes(), "dispersion")?,
        None => vec![args.dispersion; m.n_genes()],
    };

    let p = design.n_coef;
    let n = m.n_samples;
    let start = start_direction(&design);

    let per_gene = |sc: &mut Scratch, g: usize| -> Vec<f64> {
        let row = m.row(g);
        for (a, (&c, &pc)) in sc.row_aug[..n].iter_mut().zip(row.iter().zip(&prior)) {
            *a = c + pc;
        }
        let row_aug = sc.row_aug[..n].to_vec();
        let beta = fit_nb_glm(&row_aug, &design, &offset_aug, dispersions[g], &start, sc);
        beta.iter().map(|&b| b / LN2).collect()
    };

    let make = || Scratch::new(n, p);
    let rows: Vec<Vec<f64>> = if rayon::current_num_threads() > 1 {
        use rayon::prelude::*;
        (0..m.n_genes())
            .into_par_iter()
            .map_init(make, |sc, g| per_gene(sc, g))
            .collect()
    } else {
        let mut sc = make();
        (0..m.n_genes()).map(|g| per_gene(&mut sc, g)).collect()
    };

    write!(output, "{}", m.gene_col).map_err(RsomicsError::Io)?;
    for name in &design.coef_names {
        write!(output, "\t{name}").map_err(RsomicsError::Io)?;
    }
    writeln!(output).map_err(RsomicsError::Io)?;
    for (gene, coefs) in m.genes.iter().zip(&rows) {
        write!(output, "{gene}").map_err(RsomicsError::Io)?;
        for &v in coefs {
            write!(output, "\t{v:.6}").map_err(RsomicsError::Io)?;
        }
        writeln!(output).map_err(RsomicsError::Io)?;
    }
    Ok(m.n_genes() as u64)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn deviance_zero_at_mle() {
        let y = [10.0, 12.0, 8.0];
        let mu = [10.0, 12.0, 8.0];
        assert!(nb_deviance(&y, &mu, 0.1).abs() < 1e-12);
    }

    #[test]
    fn solve_diagonal() {
        let mut a = vec![2.0, 0.0, 0.0, 3.0];
        let mut b = [4.0, 9.0];
        let mut x = vec![0.0; 2];
        assert!(solve(&mut a, &mut b, &mut x, 2));
        assert!((x[0] - 2.0).abs() < 1e-12 && (x[1] - 3.0).abs() < 1e-12);
    }

    #[test]
    fn one_group_recovers_log_mean() {
        let row = [20.0, 20.0, 20.0, 20.0];
        let off = [0.0, 0.0, 0.0, 0.0];
        let b = mglm_one_group(&row, &off, 0.05);
        assert!((b - 20.0f64.ln()).abs() < 1e-9);
    }
}
