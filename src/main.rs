//! # Genetic Correlation
//!
//! This crate provides functions to compute SNP-based heritability and
//! genetic correlation between two GWAS summary statistics files.

mod cli;
mod gwas_sum;
mod mr;
use crate::cli::Args;
use crate::gwas_sum::compute_genetic_correlation;
use crate::gwas_sum::{GwasSummary, LdFile, RefAlleles};
use crate::mr::*;
use clap::Parser;
use std::io::{self};

/// Application entry point.
///
/// Parses CLI arguments and executes one of the following tasks depending on flags:
/// - vcf_clean: convert a VCF file to a summary-statistics TSV
/// - modify: modify a summary-statistics file (filter by p-value, MAF, or merge alleles)
/// - compute_h2: compute SNP heritability (requires LD file)
/// - genetic_correlation: compute genetic correlation between two sumstats (requires LD file)
/// - snp_clump: perform LD-based SNP clumping using an LD directory
/// - mr_ivm / mr_egger: run Mendelian Randomization analyses (IVW / Egger)
fn main() -> io::Result<()> {
    let args = Args::parse();

    // mission1: clean the vcf data to a summary-statistics file
    if let Some(file_path) = args.vcf_clean {
        let output_path = args.output.clone().unwrap_or_else(|| {
            // default output filename
            format!("{}.vcf_cleaned.tsv", file_path)
        });
        let gwas_vcf = GwasSummary::from_path(&file_path);
        println!(
            "cleaned the vcf file {}",
            gwas_vcf
                .filename()
                .expect("fail to read the filename of vcf file")
        );
        let snps = gwas_vcf.clean_vcf()?;
        GwasSummary::write_snps(&snps, &output_path)?;
        return Ok(());
    }

    // mission2: modify the sumstats file (maf / p / merge-allele)
    if let Some(input_path) = args.modify {
        let output_path = args.output.clone().unwrap_or_else(|| {
            // default output filename
            format!("{}.modified.tsv", input_path)
        });
        let gwas_sum = GwasSummary::from_path(&input_path);
        println!(
            "modifying the gwas summary statistics file {}",
            gwas_sum
                .filename()
                .expect("fail to read the filename of gwas sumstats")
        );
        let mut snps = gwas_sum.load_snps()?;

        if let Some(max_p) = args.max_p {
            snps = GwasSummary::filter_by_p(snps, max_p);
        }

        if let Some(snp_filter) = args.merge_alleles {
            let keep_snp_sum = RefAlleles::from_path(&snp_filter);
            let keep_snps = keep_snp_sum.load_alleles()?;
            let snp_list: Vec<(&str, &str, &str)> = RefAlleles::extract_snp_names(&keep_snps);
            let filtered_snps = GwasSummary::merge_alleles(snps, snp_list);
            snps = filtered_snps;
        }
        if let Some(min_maf) = args.min_maf {
            snps = GwasSummary::filter_by_maf(snps, min_maf);
        }

        if let Some(filter_snps_by) = args.filter_snps_by {
            snps = filter_by_snp_file(snps, &filter_snps_by);
        }

        GwasSummary::write_snps(&snps, &output_path)?;
        return Ok(());
    }

    if let Some(input_path) = args.compute_h2 {
        let gwas_sum = GwasSummary::from_path(&input_path);
        println!("computing the h2 of file {}", gwas_sum.filename().unwrap());
        if let Some(ld_path) = args.ld.clone() {
            let ld = LdFile::from_path(&ld_path);
            let h2 = gwas_sum.compute_h2(ld)?;
            println!("Heritability (h2): {}", h2)
        } else {
            panic!("please input a ld file")
        }
    }

    if let Some(paths) = args.genetic_correlation {
        if let Some(ld_path) = args.ld.clone() {
            let ld = LdFile::from_path(&ld_path);
            if paths.len() == 2 {
                let gwas_file1 = GwasSummary::from_path(&paths[0]);
                let gwas_file2 = GwasSummary::from_path(&paths[1]);
                let r_g = compute_genetic_correlation(gwas_file1, gwas_file2, ld)?;
                println!("Genetic correlation (r_g): {}", r_g);
            } else {
                eprintln!("Error: --genetic-correlation requires exactly 2 file paths");
            }
        }
    }

    // Second part: MR analysis
    if let Some(input_path) = args.snp_clump {
        let output_path = args.output.clone().unwrap_or_else(|| {
            // default output filename
            format!("{}.clumped.tsv", input_path)
        });
        let gwas_sum = GwasSummary::from_path(&input_path);
        println!("clumping snps of file {}", gwas_sum.filename().unwrap());
        if let Some(ld_dir) = args.ld_dir {
            let snps = gwas_sum.load_snps()?;
            if let Some(max_r2) = args.max_r2 {
                let clumped_snps = snp_clump(snps, max_r2, ld_dir);
                GwasSummary::write_snps(&clumped_snps, &output_path)?;
            }
        }
    }

    if args.mr_ivm || args.mr_egger {
        let exposure: GwasSummary;
        let outcome: GwasSummary;
        if let Some(exposure_path) = args.exposure {
            exposure = GwasSummary::from_path(&exposure_path);
        } else {
            panic!("please assign exposure file")
        }
        if let Some(outcome_path) = args.outcome {
            outcome = GwasSummary::from_path(&outcome_path);
        } else {
            panic!("please assign outcome file")
        }
        if args.mr_ivm {
            mr_ivm(exposure.clone(), outcome.clone())?;
        }
        if args.mr_egger {
            mr_egger(exposure, outcome)?;
        }
    }

    Ok(())
}