//! # Command Line Interface
//!
//! Defines CLI argument parsing for the GWAS utilities crate.
//! Uses `clap` to parse options for cleaning VCFs, modifying summary files,
//! computing heritability, genetic correlation, SNP clumping, and MR analyses.
use clap::Parser;

/// CLI arguments for the GWAS utilities program.
#[derive(Parser)]
#[clap(author, version, about)]
pub struct Args {
    #[clap(long)]
    pub sample_size: Option<f64>,
    /// GWAS summary statistics file to clean (VCF -> cleaned TSV)
    #[clap(long)]
    pub data_clean: Option<String>,

    /// SNP list file used to merge / filter alleles (path to allele reference file)
    #[clap(long)]
    pub merge_alleles: Option<String>,

    /// Input GWAS summary file (with Z and L2) to compute SNP-based heritability (h2)
    #[clap(long)]
    pub compute_h2: Option<String>,

    /// Path to LD scores or LD-related file required for h2 / correlation computations
    #[clap(long)]
    pub ld: Option<String>,

    /// Compute genetic correlation from two GWAS summary statistics files
    /// Provide exactly two file paths: --genetic-correlation GWAS1 GWAS2
    #[clap(long, num_args = 2, value_names = ["GWAS1", "GWAS2"])]
    pub genetic_correlation: Option<Vec<String>>,

    /// Path to a VCF file to clean (suitable for IEU OpenGWAS-style VCFs)
    #[clap(long)]
    pub vcf_clean: Option<String>,

    /// Input summary-statistics file to modify (filter by MAF / P-value / merge alleles)
    #[clap(long)]
    pub modify: Option<String>,

    /// Minimum minor allele frequency threshold for filtering (e.g. 0.01)
    #[clap(long)]
    pub min_maf: Option<f64>,

    /// Maximum p-value threshold for filtering (retain SNPs with p <= value)
    #[clap(long)]
    pub max_p: Option<f64>,

    /// Output file path for generated/modified files
    #[clap(long)]
    pub output: Option<String>,

    /// GWAS summary statistics file to perform SNP clumping on
    #[clap(long)]
    pub snp_clump: Option<String>,

    /// Directory containing per-chromosome LD r2 files (used by clumping)
    #[clap(long)]
    pub ld_dir: Option<String>,

    /// Maximum r2 threshold for SNP clumping (e.g. 0.01)
    #[clap(long)]
    pub max_r2: Option<f64>,

    /// Path to a file with SNP IDs used to filter the summary file
    #[clap(long)]
    pub filter_snps_by: Option<String>,

    /// Run IVW Mendelian Randomization (inverse-variance weighted)
    #[clap(long)]
    pub mr_ivm: bool,

    /// Run MR-Egger analysis
    #[clap(long)]
    pub mr_egger: bool,

    /// Exposure GWAS summary file for MR analyses
    #[clap(long)]
    pub exposure: Option<String>,

    /// Outcome GWAS summary file for MR analyses
    #[clap(long)]
    pub outcome: Option<String>,

    #[clap(long, num_args = 2, value_names = ["GWAS1", "GWAS2"])]
    pub pleio_identify: Option<Vec<String>>,

    #[clap(long)]
    pub min_z: Option<f64>,

    #[clap(long, num_args = 2, value_names = ["GWAS1", "GWAS2"])]
    pub coloc: Option<Vec<String>>,
}