use crate::gwas_sum::{GwasSummary, SnpInfo};
use ndarray::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{self, BufRead, BufReader};

pub fn snp_clump(snps: Vec<SnpInfo>, max_r2: f64, ld_path: String) -> Vec<SnpInfo> {
    // 1. 先按染色体分组
    let mut chr_map: HashMap<String, Vec<SnpInfo>> = HashMap::new();
    for snp in snps {
        chr_map.entry(snp.chr.clone()).or_default().push(snp);
    }

    let mut kept_snps: Vec<SnpInfo> = Vec::new();

    // 2. 对每条染色体做clump
    for (chr, chr_snps) in chr_map {
        let ld_file = format!("{}/chr{}.ld", ld_path, chr);
        let ld_pairs = match read_ld_file(&ld_file) {
            Some(map) => map,
            None => {
                // 如果没有对应LD文件，则全部保留
                kept_snps.extend(chr_snps);
                continue;
            }
        };

        // 用于标记要删除的 SNP 名称
        let mut to_remove: HashSet<String> = HashSet::new();

        // 3. 两两比对 SNP
        for i in 0..chr_snps.len() {
            for j in (i + 1)..chr_snps.len() {
                let snp1 = &chr_snps[i];
                let snp2 = &chr_snps[j];

                let key1 = (snp1.snp.clone(), snp2.snp.clone());
                let key2 = (snp2.snp.clone(), snp1.snp.clone());

                if let Some(&r2) = ld_pairs.get(&key1).or(ld_pairs.get(&key2)) {
                    if r2 > max_r2 {
                        // 比较p值，保留显著性更高的
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

        // 4. 留下未被删除的
        for snp in chr_snps.into_iter() {
            if !to_remove.contains(&snp.snp) {
                kept_snps.push(snp);
            }
        }
    }

    kept_snps
}

/// 读取LD文件：返回 HashMap<(SNP_A, SNP_B), r2>
fn read_ld_file(path: &str) -> Option<HashMap<(String, String), f64>> {
    let file = File::open(path).ok()?;
    let reader = BufReader::new(file);

    let mut map = HashMap::new();

    for (i, line) in reader.lines().enumerate() {
        let line = line.ok()?;
        if i == 0 {
            continue;
        } // 跳过表头

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

pub fn filter_by_snp_file(snps: Vec<SnpInfo>, tsv_path: &str) -> Vec<SnpInfo> {
    // 打开文件
    let file = File::open(tsv_path).expect("无法打开输入文件");
    let reader = BufReader::new(file);

    let mut snp_col_index: Option<usize> = None;
    let mut snp_set: HashSet<String> = HashSet::new();

    for (i, line) in reader.lines().enumerate() {
        let line = line.expect("读取行失败");
        let cols: Vec<&str> = line.trim().split_whitespace().collect();

        if i == 0 {
            // 自动定位 "SNP" 列
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

    // 过滤 snps
    snps.into_iter()
        .filter(|s| snp_set.contains(s.snp.as_str()))
        .collect()
}

pub fn weighted_linear_regression(
    x: &Array1<f64>,
    y: &Array1<f64>,
    w: &Array1<f64>,
) -> (f64, f64, f64, f64) {
    assert_eq!(x.len(), y.len());
    assert_eq!(x.len(), w.len());

    let w_sum = w.sum();

    // 加权均值
    let x_mean = (x * w).sum() / w_sum;
    let y_mean = (y * w).sum() / w_sum;

    // 加权协方差和方差
    let cov_xy = ((x - x_mean) * (y - y_mean) * w).sum();
    let var_x = ((x - x_mean).mapv(|v| v * v) * w).sum();

    // 斜率和截距
    let slope = cov_xy / var_x;
    let intercept = y_mean - slope * x_mean;

    // 残差
    let residuals = y - &(intercept + slope * x);

    // 加权残差平方和
    let rss = (&residuals * &residuals * w).sum();
    let df = x.len() as f64 - 2.0;
    let sigma2 = rss / df;

    // 标准误差
    let se_slope = (sigma2 / var_x).sqrt();
    let se_intercept = (sigma2 * (1.0 / w_sum + x_mean.powi(2) / var_x)).sqrt();

    (intercept, slope, se_intercept, se_slope)
}

pub fn linear_regression(x: &Array1<f64>, y: &Array1<f64>) -> (f64, f64, f64, f64) {
    assert_eq!(x.len(), y.len());
    let w = &Array1::from_elem(x.len(), 1.0);

    let w_sum = w.sum();

    // 加权均值
    let x_mean = (x * w).sum() / w_sum;
    let y_mean = (y * w).sum() / w_sum;

    // 加权协方差和方差
    let cov_xy = ((x - x_mean) * (y - y_mean) * w).sum();
    let var_x = ((x - x_mean).mapv(|v| v * v) * w).sum();

    // 斜率和截距
    let slope = cov_xy / var_x;
    let intercept = y_mean - slope * x_mean;

    // 残差
    let residuals = y - &(intercept + slope * x);

    // 加权残差平方和
    let rss = (&residuals * &residuals * w).sum();
    let df = x.len() as f64 - 2.0;
    let sigma2 = rss / df;

    // 标准误差
    let se_slope = (sigma2 / var_x).sqrt();
    let se_intercept = (sigma2 * (1.0 / w_sum + x_mean.powi(2) / var_x)).sqrt();

    (intercept, slope, se_intercept, se_slope)
}

pub fn mr_ivm(exposure: GwasSummary, outcome: GwasSummary) -> io::Result<()> {
    // 载入 SNP 列表
    let snps_exp = exposure.load_snps()?;
    let snps_out = outcome.load_snps()?; // ← 注意这里改成 outcome

    // 将 exposure SNPs 转为 HashMap，方便按 SNP 名快速匹配
    let exp_map: HashMap<String, &SnpInfo> = snps_exp
        .iter()
        .filter(|s| s.beta.is_some() && s.se.is_some())
        .map(|s| (s.snp.clone(), s))
        .collect();

    // 存储匹配后的数据
    let mut beta_exp_vec = Vec::new();
    let mut beta_out_vec = Vec::new();
    let mut weight_vec = Vec::new();

    for s_out in snps_out.iter() {
        if let (Some(beta_out), Some(se_out)) = (s_out.beta, s_out.se) {
            if let Some(s_exp) = exp_map.get(&s_out.snp) {
                if let (Some(beta_exp), Some(se_exp)) = (s_exp.beta, s_exp.se) {
                    beta_exp_vec.push(beta_exp);
                    beta_out_vec.push(beta_out);
                    // 权重为 outcome 的 1 / SE^2
                    weight_vec.push(1.0 / (se_out * se_out));
                }
            }
        }
    }

    if beta_exp_vec.is_empty() {
        return Err(io::Error::new(io::ErrorKind::Other, "没有匹配的 SNP"));
    }

    // 转换为 ndarray
    let x = Array1::from(beta_exp_vec);
    let y = Array1::from(beta_out_vec);
    let w = Array1::from(weight_vec);

    // 调用加权线性回归
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
pub fn mr_egger(exposure: GwasSummary, outcome: GwasSummary) -> io::Result<()> {
    // 载入 SNP 列表
    let snps_exp = exposure.load_snps()?;
    let snps_out = outcome.load_snps()?; // ← 注意这里改成 outcome

    // 将 exposure SNPs 转为 HashMap，方便按 SNP 名快速匹配
    let exp_map: HashMap<String, &SnpInfo> = snps_exp
        .iter()
        .filter(|s| s.beta.is_some() && s.se.is_some())
        .map(|s| (s.snp.clone(), s))
        .collect();

    // 存储匹配后的数据
    let mut beta_exp_vec = Vec::new();
    let mut beta_out_vec = Vec::new();
    let mut weight_vec = Vec::new();

    for s_out in snps_out.iter() {
        if let (Some(beta_out), Some(se_out)) = (s_out.beta, s_out.se) {
            if let Some(s_exp) = exp_map.get(&s_out.snp) {
                if let (Some(beta_exp), Some(se_exp)) = (s_exp.beta, s_exp.se) {
                    beta_exp_vec.push(beta_exp);
                    beta_out_vec.push(beta_out);
                    // 权重为 outcome 的 1 / SE^2
                    weight_vec.push(1.0 / (se_out * se_out));
                }
            }
        }
    }

    if beta_exp_vec.is_empty() {
        return Err(io::Error::new(io::ErrorKind::Other, "没有匹配的 SNP"));
    }

    // 转换为 ndarray
    let x = Array1::from(beta_exp_vec);
    let y = Array1::from(beta_out_vec);
    let w = Array1::from(weight_vec);

    // 调用加权线性回归
    let (intercept, slope, se_intercept, se_slope) = weighted_linear_regression(&x, &y, &w);

    println!("================ EGGER MR 结果 ================");
    println!("IVW 估计值(Slope)    = {:.6}", slope);
    println!("标准误(SE_slope)      = {:.6}", se_slope);
    println!("截距(Intercept)       = {:.6}", intercept);
    println!("截距标准误(SE_intercept) = {:.6}", se_intercept);
    println!("============================================");

    Ok(())
}
pub fn weighted_linear_regression_zero_intercept(
    x: &Array1<f64>,
    y: &Array1<f64>,
    w: &Array1<f64>,
) -> (f64, f64, f64, f64) {
    assert_eq!(x.len(), y.len());
    assert_eq!(x.len(), w.len());
    assert!(x.len() > 1, "Need at least 2 data points");

    let n = x.len() as f64;

    // 计算加权统计量
    let cov_xy = (x * y * w).sum();
    let var_x = (x * x * w).sum();

    assert!(var_x > 0.0, "Variance of x must be positive");

    let slope = cov_xy / var_x;

    // 残差
    let residuals = y - &(x * slope);

    // 加权残差平方和
    let rss = (&residuals * &residuals * w).sum();

    // 方法1：使用样本量-1作为自由度（更保守）
    let df = n - 1.0;
    let mse = rss / df;

    // 方法2：使用有效样本量-1（当权重差异很大时更合适）
    // let effective_n = w.sum();
    // let df = effective_n - 1.0;
    // let mse = rss / df;

    // 斜率的标准误
    let se_slope = (mse / var_x).sqrt();

    (0.0, slope, 0.0, se_slope)
}
