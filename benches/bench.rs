use std::hint::black_box;
use std::io::Write;

use criterion::{Criterion, criterion_group, criterion_main};
use rsomics_edger_predfc::{PredFcArgs, predfc};

fn make_inputs(
    n_genes: usize,
    n_samples: usize,
) -> (tempfile::TempDir, std::path::PathBuf, std::path::PathBuf) {
    let dir = tempfile::tempdir().unwrap();
    let counts = dir.path().join("counts.tsv");
    let design = dir.path().join("design.tsv");
    let mut c = std::fs::File::create(&counts).unwrap();
    write!(c, "gene").unwrap();
    for s in 0..n_samples {
        write!(c, "\tS{s}").unwrap();
    }
    writeln!(c).unwrap();
    let mut seed = 0x1234_5678u64;
    let mut rng = || {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        seed
    };
    for g in 0..n_genes {
        write!(c, "g{g}").unwrap();
        for _ in 0..n_samples {
            write!(c, "\t{}", rng() % 500).unwrap();
        }
        writeln!(c).unwrap();
    }
    let mut d = std::fs::File::create(&design).unwrap();
    writeln!(d, "Intercept\tgroup").unwrap();
    for s in 0..n_samples {
        writeln!(d, "1\t{}", if s < n_samples / 2 { 0 } else { 1 }).unwrap();
    }
    (dir, counts, design)
}

fn bench_predfc(c: &mut Criterion) {
    let (_dir, counts, design) = make_inputs(5000, 50);
    c.bench_function("predfc_5000x50", |b| {
        b.iter(|| {
            let mut out = Vec::new();
            predfc(
                &PredFcArgs {
                    counts: black_box(&counts),
                    design: black_box(&design),
                    dispersion: 0.1,
                    dispersion_file: None,
                    norm_factors: None,
                    offset_file: None,
                    prior_count: 0.125,
                },
                &mut out,
            )
            .unwrap();
            black_box(out);
        })
    });
}

criterion_group!(benches, bench_predfc);
criterion_main!(benches);
