//! Differential compat against edgeR predFC (edgeR 4.8.2, R 4.5.3).
//!
//! The committed golden was captured from `predFC(DGEList, design,
//! prior.count=0.125)`; it always runs in CI. The live-oracle test re-derives
//! the reference with Rscript when edgeR is installed and loud-skips otherwise.

use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

use rsomics_edger_predfc::{PredFcArgs, predfc};

fn golden(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

/// Parse a gene×coef TSV into (gene, values), tolerant of R's trailing-zero trim.
fn parse(text: &str) -> Vec<(String, Vec<f64>)> {
    text.lines()
        .skip(1)
        .filter(|l| !l.is_empty())
        .map(|l| {
            let mut f = l.split('\t');
            let gene = f.next().unwrap().to_string();
            let vals = f.map(|v| v.parse::<f64>().unwrap()).collect();
            (gene, vals)
        })
        .collect()
}

fn run_ours_nf() -> String {
    let mut out = Vec::new();
    predfc(
        &PredFcArgs {
            counts: &golden("counts.tsv"),
            design: &golden("design.tsv"),
            dispersion: 0.05,
            dispersion_file: Some(&golden("disp.tsv")),
            norm_factors: Some(&golden("nf.tsv")),
            offset_file: None,
            prior_count: 0.125,
        },
        &mut out,
    )
    .unwrap();
    String::from_utf8(out).unwrap()
}

fn assert_matches(reference: &str, ours: &str, tol: f64) {
    let r = parse(reference);
    let o = parse(ours);
    assert_eq!(r.len(), o.len(), "row count differs");
    for ((rg, rv), (og, ov)) in r.iter().zip(&o) {
        assert_eq!(rg, og, "gene order differs");
        assert_eq!(rv.len(), ov.len(), "coef count differs for {rg}");
        for (a, b) in rv.iter().zip(ov) {
            assert!((a - b).abs() <= tol, "gene {rg}: {a} vs {b} exceeds {tol}");
        }
    }
}

#[test]
fn matches_committed_golden() {
    let reference = std::fs::read_to_string(golden("predfc_ref.tsv")).unwrap();
    assert_matches(&reference, &run_ours_nf(), 1e-6);
}

/// A fully zero sample column drives its library size to 0, so the offset is
/// log(0) = -Inf. edgeR (edgeR 4.4.0) bails with "offsets must be finite
/// values"; we must fail loud too rather than emit an all-NaN table with a
/// success exit. The well-conditioned golden shares this design (4 samples).
#[test]
fn zero_library_column_errors() {
    let mut out = Vec::new();
    let err = predfc(
        &PredFcArgs {
            counts: &golden("counts_zerolib.tsv"),
            design: &golden("design.tsv"),
            dispersion: 0.05,
            dispersion_file: None,
            norm_factors: None,
            offset_file: None,
            prior_count: 0.125,
        },
        &mut out,
    )
    .expect_err("zero-library column must error, not emit output");
    let msg = err.to_string();
    assert!(
        msg.contains("s4") && msg.contains("offsets must be finite"),
        "unexpected error message: {msg}"
    );
    assert!(out.is_empty(), "no numeric output on the degenerate path");
}

#[test]
fn matches_live_edger() {
    let rscript = "Rscript";
    let probe = Command::new(rscript)
        .args(["-e", "stopifnot(requireNamespace('edgeR', quietly=TRUE))"])
        .status();
    match probe {
        Ok(s) if s.success() => {}
        _ => {
            eprintln!("SKIP matches_live_edger: Rscript/edgeR unavailable");
            return;
        }
    }

    let dir = std::env::temp_dir().join("rsomics-edger-predfc-compat");
    std::fs::create_dir_all(&dir).unwrap();
    let ref_path = dir.join("ref.tsv");

    let script = format!(
        r#"suppressMessages(library(edgeR))
counts <- as.matrix(read.table("{counts}", header=TRUE, row.names=1, sep="\t", check.names=FALSE))
design <- as.matrix(read.table("{design}", header=TRUE, sep="\t", check.names=FALSE))
disp <- scan("{disp}", quiet=TRUE)
nf <- scan("{nf}", quiet=TRUE)
off <- log(colSums(counts) * nf)
pf <- predFC(counts, design, prior.count=0.125, offset=off, dispersion=disp)
out <- cbind(gene=rownames(pf), as.data.frame(pf))
write.table(out, "{out}", sep="\t", quote=FALSE, row.names=FALSE)
"#,
        counts = golden("counts.tsv").display(),
        design = golden("design.tsv").display(),
        disp = golden("disp.tsv").display(),
        nf = golden("nf.tsv").display(),
        out = ref_path.display(),
    );
    let script_path = dir.join("oracle.R");
    std::fs::File::create(&script_path)
        .unwrap()
        .write_all(script.as_bytes())
        .unwrap();

    let status = Command::new(rscript).arg(&script_path).status().unwrap();
    assert!(status.success(), "oracle Rscript failed");

    let reference = std::fs::read_to_string(&ref_path).unwrap();
    assert_matches(&reference, &run_ours_nf(), 1e-5);
}
