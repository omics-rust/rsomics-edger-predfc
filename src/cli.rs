use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, Tool, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Section};

use rsomics_edger_predfc::{PredFcArgs, predfc};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-edger-predfc", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    pub counts: PathBuf,
    #[arg(long, value_name = "PATH")]
    design: PathBuf,
    #[arg(long, default_value_t = 0.05)]
    dispersion: f64,
    #[arg(long, value_name = "PATH")]
    dispersion_file: Option<PathBuf>,
    #[arg(long, value_name = "PATH")]
    norm_factors: Option<PathBuf>,
    #[arg(long, value_name = "PATH")]
    offset: Option<PathBuf>,
    #[arg(long, default_value_t = 0.125)]
    prior_count: f64,
    #[arg(short = 'o', long, default_value = "-")]
    output: String,
    #[command(flatten)]
    pub common: CommonFlags,
}

impl Tool for Cli {
    fn meta() -> ToolMeta {
        META
    }
    fn common(&self) -> &CommonFlags {
        &self.common
    }

    fn execute(self) -> Result<()> {
        self.common.install_rayon_pool()?;
        let mut out: Box<dyn std::io::Write> = if self.output == "-" {
            Box::new(std::io::stdout().lock())
        } else {
            Box::new(std::fs::File::create(&self.output).map_err(RsomicsError::Io)?)
        };
        let n = predfc(
            &PredFcArgs {
                counts: &self.counts,
                design: &self.design,
                dispersion: self.dispersion,
                dispersion_file: self.dispersion_file.as_deref(),
                norm_factors: self.norm_factors.as_deref(),
                offset_file: self.offset.as_deref(),
                prior_count: self.prior_count,
            },
            &mut out,
        )?;
        if !self.common.quiet {
            eprintln!("{n} genes");
        }
        Ok(())
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
    tagline: "Predictive (prior-count-shrunk) log2-fold-change coefficients from an NB-GLM fit (edgeR predFC).",
    origin: None,
    usage_lines: &[
        "<counts.tsv> --design design.tsv [--dispersion D | --dispersion-file f.tsv] [--norm-factors f.tsv | --offset f.tsv] [--prior-count 0.125] [-o predfc.tsv]",
    ],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[
            FlagSpec {
                short: None,
                long: "design",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("PathBuf"),
                required: true,
                default: None,
                description: "Design matrix TSV: a header of coefficient names then one numeric row per sample.",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "dispersion",
                aliases: &[],
                value: Some("<float>"),
                type_hint: Some("f64"),
                required: false,
                default: Some("0.05"),
                description: "Common negative-binomial dispersion shared across genes.",
                why_default: Some("edgeR's fallback when no per-gene dispersion is supplied."),
            },
            FlagSpec {
                short: None,
                long: "dispersion-file",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("PathBuf"),
                required: false,
                default: None,
                description: "Per-gene dispersions (one per row, gene order), overriding --dispersion.",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "norm-factors",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("PathBuf"),
                required: false,
                default: None,
                description: "Per-sample normalization factors (TMM etc.); multiplied into library sizes to form the offset.",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "offset",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("PathBuf"),
                required: false,
                default: None,
                description: "Explicit per-sample log offset, bypassing the library-size computation.",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "prior-count",
                aliases: &[],
                value: Some("<float>"),
                type_hint: Some("f64"),
                required: false,
                default: Some("0.125"),
                description: "Average prior count added (library-size scaled) before refitting.",
                why_default: Some("edgeR's default predFC shrinkage."),
            },
        ],
    }],
    examples: &[
        Example {
            description: "Predictive logFC with per-gene dispersions and TMM factors",
            command: "rsomics-edger-predfc counts.tsv --design design.tsv --dispersion-file disp.tsv --norm-factors tmm.tsv -o predfc.tsv",
        },
        Example {
            description: "Common dispersion 0.1, stronger shrinkage",
            command: "rsomics-edger-predfc counts.tsv --design design.tsv --dispersion 0.1 --prior-count 1 -o predfc.tsv",
        },
    ],
    json_result_schema_doc: None,
};

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn cli_debug_assert() {
        Cli::command().debug_assert();
    }
}
