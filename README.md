# rsomics-edger-predfc

Predictive log2-fold-change coefficients — edgeR's `predFC`. A library-size
scaled prior count is added to every count and folded into the GLM offset, the
negative-binomial GLM is refit to convergence, and its design coefficients are
returned in log2. This is the shrunken logFC edgeR recommends for ranking and
plotting, stable for genes with zero counts in a group where the raw MLE diverges.

## Usage

```
rsomics-edger-predfc counts.tsv --design design.tsv \
  --dispersion-file disp.tsv --norm-factors tmm.tsv \
  --prior-count 0.125 -o predfc.tsv
```

- `counts.tsv` — gene × sample integer matrix, first column gene id.
- `--design` — sample × coefficient design matrix, header = coefficient names.
- `--dispersion-file` — per-gene NB dispersions (e.g. from `estimateDisp`); or a
  single `--dispersion` shared across genes.
- `--norm-factors` — per-sample normalization factors (TMM etc.); the GLM offset
  is `log(lib.size · norm.factor)`. `--offset` supplies that offset directly.
- `--prior-count` — average prior count to add before refitting (default 0.125).

Output is a gene × coefficient matrix of predictive log2 fold-changes with the
design's coefficient names as the header.

## Origin

This crate is an independent Rust reimplementation of edgeR's `predFC` based on:

- McCarthy DJ, Chen Y, Smyth GK (2012). "Differential expression analysis of
  multifactor RNA-Seq experiments with respect to biological variation."
  *Nucleic Acids Research* 40(10):4288-4297. doi:10.1093/nar/gks042
- Robinson MD, McCarthy DJ, Smyth GK (2010). "edgeR: a Bioconductor package for
  differential expression analysis of digital gene expression data."
  *Bioinformatics* 26(1):139-140. doi:10.1093/bioinformatics/btp616
- The negative-binomial GLM with prior-count augmentation as documented in the
  edgeR User's Guide and `?predFC`, and black-box behavior testing against the
  edgeR binary (edgeR 4.8.2, R 4.5.3).

edgeR is GPL-licensed. No edgeR source code was read or used during
implementation — the method comes from the published papers, the package
documentation, and observed input/output behavior. Test fixtures are generated
deterministically (`set.seed`) from the edgeR pipeline.

License: MIT OR Apache-2.0.
Upstream credit: edgeR (https://bioconductor.org/packages/edgeR/), GPL (>=2).
