//! # Command Line Interface
//!
//! This module defines CLI argument parsing and main entry points
//! for running heritability, correlation, and MR analyses.
//!
//! Typically, it uses the `clap` crate to parse arguments from the command line.
use clap::Parser;

/// My genetic correlation CLI
#[derive(Parser)]
#[clap(author, version, about)]
pub struct Args {
    /// GWAS summary statistics file to clean
    #[clap(long)]
    pub data_clean: Option<String>,
    /// SNP list file for filtering
    #[clap(long)]
    pub merge_alleles: Option<String>,

    /// input the gwas summary statistics with Z and l2,compute h2 by least sqaure method.
    #[clap(long)]
    pub compute_h2: Option<String>,

    
    #[clap(long)]
    pub ld: Option<String>,

    /// Compute genetic correlation from two GWAS summary statistics
    #[clap(long, num_args = 2, value_names = ["GWAS1", "GWAS2"])]
    pub genetic_correlation: Option<Vec<String>>,

    /// clean vcf data; suitable for vcf files from IEU OpenGWAS database.
    #[clap(long)]
    pub vcf_clean: Option<String>,

    ///to filter the maf/pvalue
    #[clap(long)]
    pub modify: Option<String>,


    #[clap(long)]
    pub min_maf: Option<f64>,

    #[clap(long)]
    pub max_p: Option<f64>,

    /// Output file path
    #[clap(long)]
    pub output: Option<String>,



}
