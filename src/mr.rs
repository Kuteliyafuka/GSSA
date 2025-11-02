//! GWAS Mendelian Randomization (MR) and SNP utilities
//!
//! This module provides:
//! - SNP clumping (based on external LD files)
//! - Filtering GWAS SNPs by an external SNP list file
//! - Weighted / ordinary linear regression (including zero-intercept variant)
//! - Simple MR analysis entry points (IVW / Egger)
//!
//! Function docs explain parameters, return values and expected formats.
//!
//! Note: LD files are expected at "{ld_path}/chr{chr}.ld", with a header line.
//! Each row should contain at least 7 columns. This module uses column 3 (index 2)
//! as SNP_A, column 6 (index 5) as SNP_B, and column 7 (index 6) as r2 value.

use crate::gwas_sum::{GwasSummary, SnpInfo};
use ndarray::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{self, BufRead, BufReader};

// ...existing code...

/// Perform LD-based clumping (pruning) on a given list of SNPs.
///
/// - snps: Vector of SnpInfo to process (SnpInfo contains snp id, chr, p, beta, se, etc.)
/// - max_r2: If r2 between two SNPs > max_r2, they are considered highly linked and
///           the less significant SNP will be removed.
/// - ld_path: Directory containing LD files; function reads `{ld_path}/chr{chr}.ld`.
///
/// Returns: Vector of kept SnpInfo (original order not guaranteed).
///
/// LD file format expectation: first line is a header; each row has at least 7 columns.
/// Uses column 3 as SNP_A, column 6 as SNP_B, and column 7 as r2 (floating point).
pub fn snp_clump(snps: Vec<SnpInfo>, max_r2: f64, ld_path: String) -> Vec<SnpInfo> {
    // ...existing code...
    // 1. First, group by chromosome
    let mut chr_map: HashMap<String, Vec<SnpInfo>> = HashMap::new();
    for snp in snps {
        chr_map.entry(snp.chr.clone()).or_default().push(snp);
    }

    let mut kept_snps: Vec<SnpInfo> = Vec::new();

    // 2. Perform clumping per chromosome
    for (chr, chr_snps) in chr_map {
        let ld_file = format!("{}/chr{}.ld", ld_path, chr);
        let ld_pairs = match read_ld_file(&ld_file) {
            Some(map) => map,
            None => {
                // If no corresponding LD file, keep all SNPs
                kept_snps.extend(chr_snps);
                continue;
            }
        };

        // Used to mark SNP names to remove
        let mut to_remove: HashSet<String> = HashSet::new();

        // 3. Pairwise compare SNPs
        for i in 0..chr_snps.len() {
            for j in (i + 1)..chr_snps.len() {
                let snp1 = &chr_snps[i];
                let snp2 = &chr_snps[j];

                let key1 = (snp1.snp.clone(), snp2.snp.clone());
                let key2 = (snp2.snp.clone(), snp1.snp.clone());

                if let Some(&r2) = ld_pairs.get(&key1).or(ld_pairs.get(&key2)) {
                    if r2 > max_r2 {
                        // Compare p-values, keep the more significant one
                        let p1 = snp1.p.unwrap_or(1.0);
                        let p2 = snp2.p.unwrap_or(1.0);

                        if p1 <= p2 {
                            to_remove.insert(snp2.snp.clone());
                        } else {
                            to_remove.insert(snp1.snp.clone());
                        }
                    }
                }
            }
        }

        // 4. Keep those not marked for removal
        for snp in chr_snps.into_iter() {
            if !to_remove.contains(&snp.snp) {
                kept_snps.push(snp);
            }
        }
    }

    kept_snps
}

/// Read an LD file and return a HashMap<(SNP_A, SNP_B), r2>
///
/// - path: LD file path (string)
/// - Returns Some(map) on success, None on failure (e.g., file missing or read error).
///
/// Notes: The function skips the first line (header) and expects at least 7 columns per row.
/// If parsing r2 fails, 0.0 is used as default.
fn read_ld_file(path: &str) -> Option<HashMap<(String, String), f64>> {
    let file = File::open(path).ok()?;
    let reader = BufReader::new(file);

    let mut map = HashMap::new();

    for (i, line) in reader.lines().enumerate() {
        let line = line.ok()?;
        if i == 0 {
            continue;
        } // skip header

        let cols: Vec<&str> = line.split_whitespace().collect();
        if cols.len() < 7 {
            continue;
        }

        let snp_a = cols[2].to_string();
        let snp_b = cols[5].to_string();
        let r2: f64 = cols[6].parse().unwrap_or(0.0);

        map.insert((snp_a, snp_b), r2);
    }

    Some(map)
}

/// Filter input snps vector using SNPs listed in an external TSV/whitespace-separated file
///
/// - snps: Vector of SnpInfo to filter
/// - tsv_path: Path to a file that contains a SNP column (column name should be "SNP" or "snp")
///
/// Returns: Vector containing only SNPs that appear in the file
///
/// Note: The function prints column names when reading the header and will panic if no "SNP" column is found.
pub fn filter_by_snp_file(snps: Vec<SnpInfo>, tsv_path: &str) -> Vec<SnpInfo> {
    // Open file
    let file = File::open(tsv_path).expect("无法打开输入文件");
    let reader = BufReader::new(file);

    let mut snp_col_index: Option<usize> = None;
    let mut snp_set: HashSet<String> = HashSet::new();

    for (i, line) in reader.lines().enumerate() {
        let line = line.expect("读取行失败");
        let cols: Vec<&str> = line.trim().split_whitespace().collect();

        if i == 0 {
            // Auto-detect "SNP" column
            for (idx, col_name) in cols.iter().enumerate() {
                println!("{}", col_name);
                if col_name.eq_ignore_ascii_case("snp") {
                    snp_col_index = Some(idx);
                    break;
                }
            }
            if snp_col_index.is_none() {
                panic!("未在文件表头中找到 'SNP' 列");
            }
        } else {
            if let Some(idx) = snp_col_index {
                if let Some(value) = cols.get(idx) {
                    snp_set.insert(value.trim().to_string());
                }
            }
        }
    }

    println!("读取到 {} 个 SNP 名称。", snp_set.len());

    // Filter snps
    snps.into_iter()
        .filter(|s| snp_set.contains(s.snp.as_str()))
        .collect()
}

/// Weighted linear regression (estimates intercept and slope with weights)
///
/// - x, y: 1-D data vectors (lengths must match)
/// - w: weight vector (same length as x,y)
///
/// Returns (intercept, slope, se_intercept, se_slope):
/// - intercept: estimated intercept
/// - slope: estimated slope
/// - se_intercept: standard error of intercept
/// - se_slope: standard error of slope
///
/// The implementation uses weighted means and weighted covariance to compute estimates,
/// and estimates residual variance using residuals with degrees of freedom (n - 2).
pub fn weighted_linear_regression(
    x: &Array1<f64>,
    y: &Array1<f64>,
    w: &Array1<f64>,
) -> (f64, f64, f64, f64) {
    assert_eq!(x.len(), y.len());
    assert_eq!(x.len(), w.len());

    let w_sum = w.sum();

    // weighted mean
    let x_mean = (x * w).sum() / w_sum;
    let y_mean = (y * w).sum() / w_sum;

    // weighted covariance and variance
    let cov_xy = ((x - x_mean) * (y - y_mean) * w).sum();
    let var_x = ((x - x_mean).mapv(|v| v * v) * w).sum();

    // slope and intercept
    let slope = cov_xy / var_x;
    let intercept = y_mean - slope * x_mean;

    // residuals
    let residuals = y - &(intercept + slope * x);

    // weighted residual sum of squares
    let rss = (&residuals * &residuals * w).sum();
    let df = x.len() as f64 - 2.0;
    let sigma2 = rss / df;

    // standard errors
    let se_slope = (sigma2 / var_x).sqrt();
    let se_intercept = (sigma2 * (1.0 / w_sum + x_mean.powi(2) / var_x)).sqrt();

    (intercept, slope, se_intercept, se_slope)
}

/// Ordinary least squares linear regression (equal weights)
pub fn linear_regression(x: &Array1<f64>, y: &Array1<f64>) -> (f64, f64, f64, f64) {
    assert_eq!(x.len(), y.len());
    let w = &Array1::from_elem(x.len(), 1.0);

    let w_sum = w.sum();

    // weighted mean
    let x_mean = (x * w).sum() / w_sum;
    let y_mean = (y * w).sum() / w_sum;

    // weighted covariance and variance
    let cov_xy = ((x - x_mean) * (y - y_mean) * w).sum();
    let var_x = ((x - x_mean).mapv(|v| v * v) * w).sum();

    // slope and intercept
    let slope = cov_xy / var_x;
    let intercept = y_mean - slope * x_mean;

    // residuals
    let residuals = y - &(intercept + slope * x);

    // weighted residual sum of squares
    let rss = (&residuals * &residuals * w).sum();
    let df = x.len() as f64 - 2.0;
    let sigma2 = rss / df;

    // standard errors
    let se_slope = (sigma2 / var_x).sqrt();
    let se_intercept = (sigma2 * (1.0 / w_sum + x_mean.powi(2) / var_x)).sqrt();

    (intercept, slope, se_intercept, se_slope)
}

/// Simple IVW (inverse-variance weighted) MR analysis
///
/// - exposure / outcome: GwasSummary. The function calls load_snps() to read SNP lists and matches by SNP.
/// - Returns io::Result<()>: prints results and returns Ok(()) on success; returns Err if no matching SNPs.
///
/// Implementation details:
/// - Only uses SNPs present in both exposure and outcome with beta and se available.
/// - Weights are outcome's 1 / SE^2.
/// - Uses zero-intercept weighted linear regression (weighted_linear_regression_zero_intercept) for IVW estimate.
pub fn mr_ivm(exposure: GwasSummary, outcome: GwasSummary) -> io::Result<()> {
    // load SNP lists
    let snps_exp = exposure.load_snps()?;
    let snps_out = outcome.load_snps()?; // ← note: using outcome here

    // convert exposure SNPs to HashMap for fast matching by SNP id
    let exp_map: HashMap<String, &SnpInfo> = snps_exp
        .iter()
        .filter(|s| s.beta.is_some() && s.se.is_some())
        .map(|s| (s.snp.clone(), s))
        .collect();

    // store matched data
    let mut beta_exp_vec = Vec::new();
    let mut beta_out_vec = Vec::new();
    let mut weight_vec = Vec::new();

    for s_out in snps_out.iter() {
        if let (Some(beta_out), Some(se_out)) = (s_out.beta, s_out.se) {
            if let Some(s_exp) = exp_map.get(&s_out.snp) {
                if let (Some(beta_exp), Some(_se_exp)) = (s_exp.beta, s_exp.se) {
                    beta_exp_vec.push(beta_exp);
                    beta_out_vec.push(beta_out);
                    // weight is outcome's 1 / SE^2
                    weight_vec.push(1.0 / (se_out * se_out));
                }
            }
        }
    }

    if beta_exp_vec.is_empty() {
        return Err(io::Error::new(io::ErrorKind::Other, "没有匹配的 SNP"));
    }

    // convert to ndarray
    let x = Array1::from(beta_exp_vec);
    let y = Array1::from(beta_out_vec);
    let w = Array1::from(weight_vec);

    // call weighted linear regression
    let (intercept, slope, se_intercept, se_slope) =
        weighted_linear_regression_zero_intercept(&x, &y, &w);

    println!("================ IVM MR 结果 ================");
    println!("IVW 估计值(Slope)    = {:.6}", slope);
    println!("标准误(SE_slope)      = {:.6}", se_slope);
    println!("截距(Intercept)       = {:.6}", intercept);
    println!("截距标准误(SE_intercept) = {:.6}", se_intercept);
    println!("============================================");

    Ok(())
}

/// Simple MR-Egger analysis entry point
///
/// - Similar to mr_ivm, but uses weighted linear regression with intercept (weighted_linear_regression)
/// - Prints the estimated intercept (used to detect directional pleiotropy)
/// - Returns io::Result<()>
pub fn mr_egger(exposure: GwasSummary, outcome: GwasSummary) -> io::Result<()> {
    // load SNP lists
    let snps_exp = exposure.load_snps()?;
    let snps_out = outcome.load_snps()?; // ← note: using outcome here

    // convert exposure SNPs to HashMap for fast matching by SNP id
    let exp_map: HashMap<String, &SnpInfo> = snps_exp
        .iter()
        .filter(|s| s.beta.is_some() && s.se.is_some())
        .map(|s| (s.snp.clone(), s))
        .collect();

    // store matched data
    let mut beta_exp_vec = Vec::new();
    let mut beta_out_vec = Vec::new();
    let mut weight_vec = Vec::new();

    for s_out in snps_out.iter() {
        if let (Some(beta_out), Some(se_out)) = (s_out.beta, s_out.se) {
            if let Some(s_exp) = exp_map.get(&s_out.snp) {
                if let (Some(beta_exp), Some(_se_exp)) = (s_exp.beta, s_exp.se) {
                    beta_exp_vec.push(beta_exp);
                    beta_out_vec.push(beta_out);
                    // weight is outcome's 1 / SE^2
                    weight_vec.push(1.0 / (se_out * se_out));
                }
            }
        }
    }

    if beta_exp_vec.is_empty() {
        return Err(io::Error::new(io::ErrorKind::Other, "没有匹配的 SNP"));
    }

    // convert to ndarray
    let x = Array1::from(beta_exp_vec);
    let y = Array1::from(beta_out_vec);
    let w = Array1::from(weight_vec);

    // call weighted linear regression
    let (intercept, slope, se_intercept, se_slope) = weighted_linear_regression(&x, &y, &w);

    println!("================ EGGER MR 结果 ================");
    println!("IVW 估计值(Slope)    = {:.6}", slope);
    println!("标准误(SE_slope)      = {:.6}", se_slope);
    println!("截距(Intercept)       = {:.6}", intercept);
    println!("截距标准误(SE_intercept) = {:.6}", se_intercept);
    println!("============================================");

    Ok(())
}

/// Zero-intercept (through-origin) weighted linear regression implementation
///
/// - Returns (intercept, slope, se_intercept, se_slope);
///   intercept and se_intercept are set to 0.0 in the zero-intercept model
///   to keep the same function signature.
///
/// Details: uses weighted <x*y> / <x^2> as slope estimate; residual variance
/// is estimated using n-1 degrees of freedom (more conservative).
pub fn weighted_linear_regression_zero_intercept(
    x: &Array1<f64>,
    y: &Array1<f64>,
    w: &Array1<f64>,
) -> (f64, f64, f64, f64) {
    assert_eq!(x.len(), y.len());
    assert_eq!(x.len(), w.len());
    assert!(x.len() > 1, "Need at least 2 data points");

    let n = x.len() as f64;

    // compute weighted statistics
    let cov_xy = (x * y * w).sum();
    let var_x = (x * x * w).sum();

    assert!(var_x > 0.0, "Variance of x must be positive");

    let slope = cov_xy / var_x;

    // residuals
    let residuals = y - &(x * slope);

    // weighted residual sum of squares
    let rss = (&residuals * &residuals * w).sum();

    // Method 1: use n-1 degrees of freedom (more conservative)
    let df = n - 1.0;
    let mse = rss / df;

    // Method 2: use effective sample size - 1 (better when weights vary greatly)
    // let effective_n = w.sum();
    // let df = effective_n - 1.0;
    // let mse = rss / df;

    // standard error of slope
    let se_slope = (mse / var_x).sqrt();

    (0.0, slope, 0.0, se_slope)
}