//! GWAS summary utilities and LD file helpers
//!
//! This module provides:
//! - GwasSummary: load/clean GWAS summary-statistics files and basic filters
//! - SnpInfo / AllelesInfo: per-SNP data structures and parsers
//! - LdFile: load LD score mappings
//! - compute_genetic_correlation: estimate genetic correlation between two GWAS files
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use ndarray::prelude::*;
// ...existing code...

/// Represents a GWAS summary-statistics file.
#[derive(Debug, Clone)]
pub struct GwasSummary {
    /// Path to the summary-statistics file.
    pub summary_path: String,
}

/// Represents a reference alleles file.
#[derive(Debug, Clone)]
pub struct RefAlleles {
    /// Path to the alleles file.
    pub alleles_path: String,
}

/// Stores allele information for a single SNP.
#[derive(Debug, Clone)]
pub struct AllelesInfo {
    /// SNP identifier (e.g. "rs123456").
    pub snp: String,
    /// Allele 1 (e.g. "A").
    pub a1: String,
    /// Allele 2 (e.g. "G").
    pub a2: String,
}

#[derive(Debug, Clone)]
pub struct LdFile {
    /// Path to the LD score file.
    pub ld_path: String,
}

/// Stores information for one SNP in GWAS summary statistics.
#[derive(Debug, Clone)]
pub struct SnpInfo {
    /// SNP identifier (e.g. "rs123456")
    pub snp: String,
    /// Chromosome name/number
    pub chr: String,
    /// Position (bp)
    pub pos: u32,
    /// Alleles (A1 / A2)
    pub a1: String,
    pub a2: String,
    /// Minor allele frequency (optional)
    pub maf: Option<f64>,
    /// Effect size (beta or log-OR)
    pub beta: Option<f64>,
    /// Standard error of effect size
    pub se: Option<f64>,
    /// P-value
    pub p: Option<f64>,
    /// Effective sample size
    pub neff: Option<f64>,
}

impl SnpInfo {
    /// Parse a TSV line produced by `clean_vcf()` into SnpInfo.
    ///
    /// Expected columns (tab-separated):
    /// CHR, POS, SNP, A1, A2, BETA, SE, P, MAF, NEFF
    pub fn from_cleaned_line(line: &str) -> Option<Self> {
        if line.starts_with('#') || line.trim().is_empty() {
            return None;
        }

        let fields: Vec<&str> = line.trim().split('\t').collect();

        if fields.len() < 10 {
            panic!("Input file format error: not enough columns. Ensure the file was produced by clean_vcf.");
        }

        Some(Self {
            chr: fields[0].to_string(),
            pos: fields[1].parse().ok()?,
            snp: fields[2].to_string(),
            a1: fields[3].to_string(),
            a2: fields[4].to_string(),
            beta: fields[5].parse().ok(),
            se: fields[6].parse().ok(),
            p: fields[7].parse().ok(),
            maf: fields[8].parse().ok(),
            neff: fields[9].parse().ok(),
        })
    }

    /// Parse a VCF-like line into SnpInfo.
    ///
    /// The function expects at least 10 columns and that the FORMAT/SAMPLE fields
    /// include tags such as ES, SE, LP, AF, SS when present.
    pub fn from_vcf_line(line: &str) -> Option<Self> {
        if line.starts_with('#') {
            return None;
        }

        let fields: Vec<&str> = line.trim().split('\t').collect();
        if fields.len() < 10 {
            return None;
        }

        let chr = fields[0].to_string();
        let pos = fields[1].parse::<u32>().ok()?;
        let snp = fields[2].to_string();
        let a1 = fields[3].to_string();
        let a2 = fields[4].to_string();

        // Parse FORMAT and sample fields (colon-separated key:value)
        let format = fields[8];
        let values = fields[9];
        let keys: Vec<&str> = format.split(':').collect();
        let vals: Vec<&str> = values.split(':').collect();
        if keys.len() != vals.len() {
            return None;
        }

        let map: HashMap<_, _> = keys.iter().zip(vals.iter()).collect();

        let beta = map.get(&"ES").and_then(|v| v.parse::<f64>().ok());
        let se = map.get(&"SE").and_then(|v| v.parse::<f64>().ok());
        let lp = map.get(&"LP").and_then(|v| v.parse::<f64>().ok());
        let maf = map.get(&"AF").and_then(|v| v.parse::<f64>().ok());
        let neff = map.get(&"SS").and_then(|v| v.parse::<f64>().ok());

        let p = lp.map(|v| 10f64.powf(-v));

        Some(Self {
            snp,
            chr,
            pos,
            a1,
            a2,
            maf,
            beta,
            se,
            p,
            neff,
        })
    }

    /// Format SnpInfo as a TSV line matching the header written by `write_snps`.
    pub fn to_tsv(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.chr,
            self.pos,
            self.snp,
            self.a1,
            self.a2,
            self.beta.map_or("NA".to_string(), |v| v.to_string()),
            self.se.map_or("NA".to_string(), |v| v.to_string()),
            self.p.map_or("NA".to_string(), |v| format!("{:.6e}", v)),
            self.maf.map_or("NA".to_string(), |v| v.to_string()),
            self.neff.map_or("NA".to_string(), |v| v.to_string()),
        )
    }
}

impl GwasSummary {
    /// Create a GwasSummary that points to `path`.
    pub fn from_path(path: &str) -> Self {
        Self {
            summary_path: path.to_string(),
        }
    }

    /// Read a VCF-like file and extract cleaned SNP records.
    ///
    /// Returns a vector of SnpInfo parsed from sample columns.
    pub fn clean_vcf(&self) -> io::Result<Vec<SnpInfo>> {
        let file = File::open(&self.summary_path)?;
        let reader = BufReader::new(file);
        let mut snps = Vec::new();
        for line in reader.lines() {
            let line = line?;
            if let Some(snp) = SnpInfo::from_vcf_line(&line) {
                snps.push(snp);
            }
        }
        Ok(snps)
    }

    /// Load SNPs from a cleaned GWAS summary file (TSV).
    ///
    /// Skips empty lines and the header line that contains "CHR".
    pub fn load_snps(&self) -> io::Result<Vec<SnpInfo>> {
        let file = File::open(&self.summary_path)?;
        let reader = BufReader::new(file);
        let mut snps = Vec::new();
        for (_i, line) in reader.lines().enumerate() {
            let line = line?;
            // Skip empty lines or header
            if line.trim().is_empty() {
                continue;
            }
            if line.contains("CHR") {
                continue; // skip header row
            }

            if let Some(snp_info) = SnpInfo::from_cleaned_line(&line) {
                snps.push(snp_info);
            }
        }

        Ok(snps)
    }

    /// Return the filename component of the summary_path, if available.
    pub fn filename(&self) -> Option<String> {
        std::path::Path::new(&self.summary_path)
            .file_name()
            .and_then(|s| s.to_str())
            .map(|s| s.to_string())
    }

    /// Filter SNPs by p-value threshold (retain SNPs with p <= threshold).
    pub fn filter_by_p(snps: Vec<SnpInfo>, p_threshold: f64) -> Vec<SnpInfo> {
        let result: Vec<SnpInfo> = snps
            .into_iter()
            .filter(|s| s.p.map_or(false, |p| p <= p_threshold))
            .collect();
        println!(
            "Filtered by p <= {}; remaining SNPs: {}",
            p_threshold,
            result.len()
        );
        result
    }


    /// Filter SNPs by MAF threshold (retain SNPs with maf >= threshold).
    pub fn filter_by_maf(snps: Vec<SnpInfo>, maf_threshold: f64) -> Vec<SnpInfo> {
        let result: Vec<SnpInfo> = snps
            .into_iter()
            .filter(|s| s.maf.map_or(false, |maf| maf >= maf_threshold && maf <= 1.0 -maf_threshold))
            .collect();
        println!(
            "Filtered by MAF >= {}; remaining SNPs: {}",
            maf_threshold,
            result.len()
        );
        result
    }

    pub fn add_sample_size(snps: Vec<SnpInfo>, sample_size:f64) -> Vec<SnpInfo> {
        let result: Vec<SnpInfo> = snps
            .into_iter()
            .map(|mut s| {
                s.neff = Some(sample_size); 
                s})
            .collect();
        println!(
            "set N = {}",
            sample_size
        );
        result
    }

    /// Keep SNPs whose (SNP, A2, A1) or (SNP, A1, A2) appear in `snp_list`.
    ///
    /// When an entry matches with swapped alleles, flip MAF (maf -> 1 - maf).
    pub fn merge_alleles<'a>(
        snps: Vec<SnpInfo>,
        snp_list: Vec<(&'a str, &'a str, &'a str)>,
    ) -> Vec<SnpInfo> {
        let mut flip_count = 0;
        let set: HashSet<(&str, &str, &str)> = snp_list.into_iter().collect();
        let mut result1: Vec<SnpInfo> = snps
            .clone()
            .into_iter()
            .filter(|s| {
                let key = (s.snp.as_str(), s.a2.as_str(), s.a1.as_str());
                set.contains(&key)
            })
            .map(|mut s| {
                if let Some(beta) = s.beta {
                    s.beta = Some(0.0 - beta);
                }
                let a1 = s.a1;
                s.a1 = s.a2;
                s.a2 = a1;
                flip_count += 1;
                s
            })
            .collect();
        let result2: Vec<SnpInfo> = snps
            .into_iter()
            .filter(|s| {
                let key = (s.snp.as_str(), s.a1.as_str(), s.a2.as_str());
                set.contains(&key)
            })
            .map(|mut s| {
                if let Some(maf) = s.maf {
                    s.maf = Some(1.0 - maf);
                }
                s
            })
            .collect();
        result1.extend(result2);
        println!("Merged alleles; remaining SNPs: {},flipped alleles:{}", result1.len(),flip_count);
        result1
    }

    /// Compute SNP-based heritability (h^2) using LD scores.
    ///
    /// Requires LD scores mapped by SNP name via `LdFile`.
    pub fn compute_h2(&self, ld: LdFile) -> io::Result<f64> {
        // 1. Load SNPs from the GWAS summary
        let snps = self.load_snps()?; // Vec<SnpInfo>

        // 2. Load LD scores: HashMap<SNP, ldscore>
        let ld_map = ld.load_ld_scores()?; // HashMap<String, f64>

        // 3. Build vectors for regression: x = L2, y = Z^2 - 1, collect neff
        let mut x_vec = Vec::new();
        let mut y_vec = Vec::new();
        let mut neff_vec = Vec::new();
        let mut w_vec = Vec::new();
        for snp in snps.iter() {
            let beta = match snp.beta { Some(b) => b, None => continue };
            let se = match snp.se { Some(b) => b, None => continue };
            let neff = match snp.neff { Some(b) => b, None => continue };
            let z = beta / se;

            let l2 = match ld_map.get(&snp.snp) {
                Some(&val) => val,
                None => continue,
            };

            x_vec.push(l2);
            y_vec.push(z * z - 1.0);
            w_vec.push(1.0 / (se * se));  // ✅ 保证与 x_vec/y_vec 同步
            neff_vec.push(neff);
        }
        let m = x_vec.len() as f64; // number of SNPs
        let mean_neff = neff_vec.iter().sum::<f64>() / m;
        println!("x_vec len = {}, y_vec len = {}, neff_vec len = {}", 
            x_vec.len(), y_vec.len(), neff_vec.len());
        // // Ordinary least squares slope
        // let x_mean = x_vec.iter().sum::<f64>() / m;
        // let y_mean = y_vec.iter().sum::<f64>() / m;

        // let numerator: f64 = x_vec
        //     .iter()
        //     .zip(y_vec.iter())
        //     .map(|(x, y)| (x - x_mean) * (y - y_mean))
        //     .sum();

        // let denominator: f64 = x_vec.iter().map(|x| (x - x_mean).powi(2)).sum();
        // let slope = numerator / denominator;

        // h^2 = slope * M / mean_neff

        /*let x_arr: Array1<f64> = Array1::from(x_vec);
        let y_arr: Array1<f64> = Array1::from(y_vec);
        let (intercept, slope, se_intercept, se_slope) = linear_regression(&x_arr, &y_arr);
        let h2 = slope * m / mean_neff;
        println!(
            "h2 {:.4},slope {:.4},intercept {:.4}, SNPs {}, mean effective sample size {},se_slope {:.4}, se_intercept {:.4}",
            h2, slope, intercept, m, mean_neff, se_slope, se_intercept
        );
        Ok(h2)*/
        // 构造 ndarray
        let x_arr: Array1<f64> = Array1::from(x_vec);
        let y_arr: Array1<f64> = Array1::from(y_vec);

        // 权重：用每个 SNP 的 se 计算权重 w = 1 / se^2
        //let w_vec: Vec<f64> = neff_vec.iter().map(|&n| n).collect(); // 或者换成 1/se^2
        let w_arr: Array1<f64> = Array1::from(w_vec);

        // 加权线性回归
        let w_sum = w_arr.sum();
        let x_mean = (&x_arr * &w_arr).sum() / w_sum;
        let y_mean = (&y_arr * &w_arr).sum() / w_sum;

        let numerator = ((&x_arr - x_mean) * (&y_arr - y_mean) * &w_arr).sum();
        let denominator = ((&x_arr - x_mean).mapv(|v| v * v) * &w_arr).sum();

        let slope = numerator / denominator;
        let intercept = y_mean - slope * x_mean;

        // 残差、加权标准误计算
        let residuals = &y_arr - &(intercept + &(slope * &x_arr));
        let sigma2 = (&w_arr * &(&residuals * &residuals)).sum() / (w_sum - 2.0);
        let se_slope = (sigma2 / denominator).sqrt();
        let se_intercept = (sigma2 * (1.0 / w_sum + x_mean.powi(2) / denominator)).sqrt();

        // 估算 h²
        let h2 = slope * m / mean_neff;

        println!(
            "Weighted h2 {:.4}, slope {:.4}, intercept {:.4}, SNPs {}, mean effective sample size {}, se_slope {:.4}, se_intercept {:.4}",
            h2, slope, intercept, m, mean_neff, se_slope, se_intercept
        );

        Ok(h2)
    }

    /// Write a list of SnpInfo to a TSV file (header included).
    pub fn write_snps(snps: &[SnpInfo], output_path: &str) -> io::Result<()> {
        let output_file = File::create(output_path)?;
        let mut writer = BufWriter::new(output_file);
        writeln!(writer, "CHR\tPOS\tSNP\tA1\tA2\tBETA\tSE\tP\tMAF\tN")?;
        for snp in snps {
            writeln!(writer, "{}", snp.to_tsv())?;
        }
        Ok(())
    }
}

impl LdFile {
    /// Create an LdFile pointing to `path`.
    pub fn from_path(path: &str) -> Self {
        Self {
            ld_path: path.to_string(),
        }
    }

    /// Load LD scores from a file into a HashMap<SNP, ld_score>.
    ///
    /// The function skips a header line containing "SNP" and expects the LD score
    /// to be in column index 5 (sixth column).
    pub fn load_ld_scores(&self) -> io::Result<HashMap<String, f64>> {
        let file = File::open(&self.ld_path)?;
        let reader = BufReader::new(file);
        let mut ld_map = HashMap::new();
        // skip header
        for (_i, line) in reader.lines().enumerate() {
            let line = line?;
            if line.contains("SNP") {
                continue;
            }
            let fields: Vec<&str> = line.split_whitespace().collect();
            if fields.len() < 6 {
                continue; // skip incomplete lines
            }
            let snp = fields[1].to_string(); // SNP column
            let ld_score = fields[5].parse::<f64>().unwrap_or(0.0); // L2 column
            ld_map.insert(snp, ld_score);
        }

        Ok(ld_map)
    }
}

impl RefAlleles {
    /// Create a RefAlleles that points to `path`.
    pub fn from_path(path: &str) -> Self {
        Self {
            alleles_path: path.to_string(),
        }
    }

    /// Load allele reference file into a vector of AllelesInfo.
    ///
    /// Skips empty lines and header lines containing "CHR".
    pub fn load_alleles(&self) -> io::Result<Vec<AllelesInfo>> {
        let file = File::open(&self.alleles_path)?;
        let reader = BufReader::new(file);
        let mut snps = Vec::new();
        for (_i, line) in reader.lines().enumerate() {
            let line = line?;
            // skip empty lines or header
            if line.trim().is_empty() {
                continue;
            }
            if line.contains("CHR") {
                continue; // skip header row
            }
            if let Some(snp_info) = AllelesInfo::from_cleaned_line(&line) {
                snps.push(snp_info);
            }
        }
        Ok(snps)
    }

    /// Extract (SNP, A2, A1) tuples from AllelesInfo slice for matching.
    pub fn extract_snp_names(snps: &[AllelesInfo]) -> Vec<(&str, &str, &str)> {
        snps.iter()
            .map(|s| (s.snp.as_str(), s.a1.as_str(), s.a2.as_str()))
            .collect()
    }
}

impl AllelesInfo {
    /// Parse a cleaned allele line into AllelesInfo.
    ///
    /// Expected columns (tab-separated): SNP, A1, A2
    pub fn from_cleaned_line(line: &str) -> Option<Self> {
        if line.starts_with('#') || line.trim().is_empty() {
            return None;
        }
        let fields: Vec<&str> = line.trim().split('\t').collect();

        if fields.len() < 3 {
            panic!("Input alleles file format error: not enough columns. Ensure the file was produced by clean_vcf.");
        }

        Some(Self {
            snp: fields[0].to_string(),
            a1: fields[1].to_string(),
            a2: fields[2].to_string(),
        })
    }
}

/// Compute genetic correlation between two GWAS summary files using LD scores.
///
/// Returns an estimate of r_g (genetic correlation).
pub fn compute_genetic_correlation(
    file_path1: GwasSummary,
    file_path2: GwasSummary,
    ld_path: LdFile,
) -> io::Result<f64> {
    // compute SNP heritability for both traits
    let h1_2 = file_path1.compute_h2(ld_path.clone())?;
    let h2_2 = file_path2.compute_h2(ld_path.clone())?;
    let ld_map = ld_path.load_ld_scores()?;
    let h1_h2 = (h1_2 * h2_2).sqrt();
    let mut map1_nol: HashMap<String, (f64, f64)> = HashMap::new();
    let mut map2_nol: HashMap<String, (f64, f64)> = HashMap::new();

    let snps1 = file_path1.load_snps()?;
    let snps2 = file_path2.load_snps()?;
    for snp in snps1.iter() {
        // require beta, se and neff
        let beta = match snp.beta {
            Some(b) => b,
            None => continue,
        };
        let se = match snp.se {
            Some(b) => b,
            None => continue,
        };
        let neff = match snp.neff {
            Some(b) => b,
            None => continue,
        };
        let z = beta / se;

        map1_nol.insert(snp.snp.clone(), (z, neff));
    }

    for snp in snps2.iter() {
        // require beta, se and neff
        let beta = match snp.beta {
            Some(b) => b,
            None => continue,
        };
        let se = match snp.se {
            Some(b) => b,
            None => continue,
        };
        let neff = match snp.neff {
            Some(b) => b,
            None => continue,
        };
        let z = beta / se;

        map2_nol.insert(snp.snp.clone(), (z, neff));
    }

    // attach LD scores to z / neff entries, dropping SNPs without LD score
    let map1: HashMap<String, (f64, f64, f64)> = map1_nol
        .into_iter()
        .flat_map(|(id, (z, neff))| ld_map.get(&id).map(|&l2| (id, (z, l2, neff))))
        .collect();
    let map2: HashMap<String, (f64, f64, f64)> = map2_nol
        .into_iter()
        .flat_map(|(id, (z, neff))| ld_map.get(&id).map(|&l2| (id, (z, l2, neff))))
        .collect();

    let keys1: HashSet<_> = map1.keys().cloned().collect();
    let keys2: HashSet<_> = map2.keys().cloned().collect();
    let inter_keys: HashSet<_> = keys1.intersection(&keys2).cloned().collect();
    let map1_filtered: HashMap<String, (f64, f64, f64)> = map1
        .into_iter()
        .filter(|(k, _)| inter_keys.contains(k))
        .collect();
    let map2_filtered: HashMap<String, (f64, f64, f64)> = map2
        .into_iter()
        .filter(|(k, _)| inter_keys.contains(k))
        .collect();

    let n1_vec: Vec<f64> = map1_filtered.values().map(|(_, _, neff)| *neff).collect();
    let n2_vec: Vec<f64> = map2_filtered.values().map(|(_, _, neff)| *neff).collect();
    let m = n1_vec.len() as f64;
    let n1 = n1_vec.iter().sum::<f64>() / m;
    let n2 = n2_vec.iter().sum::<f64>() / m;

    // build vectors z1*z2 and L2
    let mut z1z2_vec: Vec<f64> = Vec::new();
    let mut l2_vec: Vec<f64> = Vec::new();

    for key in &inter_keys {
        let (z1, l2_1, _) = map1_filtered.get(key).unwrap();
        let (z2, l2_2, _) = map2_filtered.get(key).unwrap();
        z1z2_vec.push(z1 * z2);
        // take average L2 for robustness
        l2_vec.push((l2_1 + l2_2) / 2.0);
    }

    // // OLS slope for regression of z1*z2 on L2
    // let x_mean = l2_vec.iter().sum::<f64>() / m;
    // let y_mean = z1z2_vec.iter().sum::<f64>() / m;

    // let numerator: f64 = l2_vec
    //     .iter()
    //     .zip(z1z2_vec.iter())
    //     .map(|(x, y)| (x - x_mean) * (y - y_mean))
    //     .sum();

    // let denominator: f64 = l2_vec.iter().map(|x| (x - x_mean).powi(2)).sum();

    // let slope = numerator / denominator;

    // theoretical formula to convert slope to genetic correlation
    /*let x_arr: Array1<f64> = Array1::from(l2_vec);
    let y_arr: Array1<f64> = Array1::from(z1z2_vec);
    let (intercept, slope, _se_intercept, _se_slope) = linear_regression(&x_arr, &y_arr);
    let genetic_correlation = slope * m / ((n1 * n2).sqrt() * h1_h2);

    println!(
        "M = {:.0}, N1 = {:.0}, N2 = {:.0}, h1² = {:.4}, h2² = {:.4}, slope = {:.6}, r_g = {:.6},intercept = {:.4}",
        m, n1, n2, h1_2, h2_2, slope, genetic_correlation, intercept
    );*/
    // 计算权重
    let weights: Vec<f64> = l2_vec.iter().map(|l| 1.0 / (l + 1.0).powi(2)).collect();

    // 加权平均
    let w_sum: f64 = weights.iter().sum();
    let x_mean = l2_vec.iter().zip(weights.iter()).map(|(x, w)| x * w).sum::<f64>() / w_sum;
    let y_mean = z1z2_vec.iter().zip(weights.iter()).map(|(y, w)| y * w).sum::<f64>() / w_sum;

    // 加权协方差与方差
    let cov_xy = l2_vec.iter()
        .zip(z1z2_vec.iter())
        .zip(weights.iter())
        .map(|((x, y), w)| w * (x - x_mean) * (y - y_mean))
        .sum::<f64>();

    let var_x = l2_vec.iter()
        .zip(weights.iter())
        .map(|(x, w)| w * (x - x_mean).powi(2))
        .sum::<f64>();

    let slope = cov_xy / var_x;
    let intercept = y_mean - slope * x_mean;

    // 计算遗传相关
    let rg = slope * m / ((n1 * n2).sqrt() * h1_h2);
    println!(
        "M = {:.0}, N1 = {:.0}, N2 = {:.0}, h1² = {:.4}, h2² = {:.4}, slope = {:.6}, r_g = {:.6},intercept = {:.4}",
        m, n1, n2, h1_2, h2_2, slope, rg, intercept
    );
    Ok(rg)
}