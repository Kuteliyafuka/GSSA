//! # Genetic Correlation
//!
//! This crate provides functions to compute SNP-based heritability and
//! genetic correlation between two GWAS summary statistics files.

mod cli;
mod gwas_sum;
use crate::cli::Args;
use crate::gwas_sum::{GwasSummary, LdFile, RefAlleles};
use clap::Parser;
use std::io::{self};
fn main() -> io::Result<()> {
    let args = Args::parse();

    /*if let Some(paths) = args.genetic_correlation {
        if paths.len() == 2 {
            let file1 = &paths[0];
            let file2 = &paths[1];
            let r_g = compute_genetic_correlation(file1, file2)?;
            println!("Genetic correlation (r_g): {}", r_g);
        } else {
            eprintln!("Error: --genetic-correlation requires exactly 2 file paths");
        }
    }*/

    //mission1: clean the vcf data to sumstats file
    if let Some(file_path) = args.vcf_clean {
        let output_path = args.output.clone().unwrap_or_else(|| {
            // 默认输出文件名
            format!("{}.vcf_cleaned.tsv", file_path)
        });
        let gwas_vcf = GwasSummary::from_path(&file_path);
        println!("cleaned the vcf file {}",gwas_vcf.filename().expect("fail to read the filename of vcf file"));
        let snps = gwas_vcf.clean_vcf()?;
        GwasSummary::write_snps(&snps, &output_path)?;
        return Ok(());
    }


    //mission2: modify the sumstats file maf/p/merge-allele
    if let Some(input_path) = args.modify {
        let output_path = args.output.clone().unwrap_or_else(|| {
            // 默认输出文件名
            format!("{}.modified.tsv", input_path)
        });
        let gwas_sum = GwasSummary::from_path(&input_path);
        println!("modifyed the gwas summary statistics file {}",gwas_sum.filename().expect("fail to read the filename of gwas sumstats"));
        let mut snps = gwas_sum.load_snps()?;

        if let Some(max_p) = args.max_p {
            snps = GwasSummary::filter_by_p(snps, max_p);
        }
        if let Some(min_maf) = args.min_maf {
            snps = GwasSummary::filter_by_maf(snps, min_maf);
        }

        if let Some(snp_filter) = args.merge_alleles {
            let keep_snp_sum =  RefAlleles::from_path(&snp_filter);
            let keep_snps = keep_snp_sum.load_alleles()?;
            let snp_list:Vec<(String,String,String)> = RefAlleles::extract_snp_names(&keep_snps);
            let filtered_snps = GwasSummary::merge_alleles(snps, snp_list);
            snps = filtered_snps;
        }
        GwasSummary::write_snps(&snps, &output_path)?;
        return Ok(());
    }
    

    if let Some(input_path) = args.compute_h2 {
        let gwas_sum = GwasSummary::from_path(&input_path);
        println!("computed the h2 of file {}",gwas_sum.filename().unwrap());
        if let Some(ld_path) = args.ld {
            let ld = LdFile::from_path(&ld_path);
            let h2 = gwas_sum.compute_h2(ld)?;
            println!("遗传率是{}",h2)
        } else {
            panic!("please input a ld file")
        }
    }

    Ok(())
}