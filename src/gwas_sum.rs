use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};

/// 存储GWAS汇总统计信息文件的信息
#[derive(Debug, Clone)]
pub struct GwasSummary {
    /// 汇总统计文件路径
    pub summary_path: String,
}
/// 存储GWAS汇总统计信息文件的信息
#[derive(Debug, Clone)]
pub struct RefAlleles {
    /// 汇总统计文件路径
    pub alleles_path: String,
}

/// 存储单个SNP（单核苷酸多态性）的信息
#[derive(Debug, Clone)]
pub struct AllelesInfo {
    /// SNP名称（如 rs123456）
    pub snp: String,

    /// 等位基因（如 "A" / "G"）
    pub a1: String,
    pub a2: String,
}

#[derive(Debug, Clone)]
pub struct LdFile {
    /// 汇总统计文件路径
    pub ld_path: String,
}
/// 存储单个SNP（单核苷酸多态性）的信息
#[derive(Debug, Clone)]
pub struct SnpInfo {
    /// SNP名称（如 rs123456）
    pub snp: String,
    /// 所在染色体编号
    pub chr: String,
    /// 在染色体上的位置（bp）
    pub pos: u32,
    /// 等位基因（如 "A" / "G"）
    pub a1: String,
    pub a2: String,
    /// 频率（如次要等位基因频率）
    pub maf: Option<f64>,
    /// 效应值（beta 或 OR）
    pub beta: Option<f64>,
    /// 标准误
    pub se: Option<f64>,
    /// P 值
    pub p: Option<f64>,
    /// effective sample size
    pub neff: Option<f64>,
}

impl SnpInfo {
    pub fn from_cleaned_line(line: &str) -> Option<Self> {
        if line.starts_with('#') || line.trim().is_empty() {
            return None;
        }

        let fields: Vec<&str> = line.trim().split('\t').collect();

        if fields.len() < 10 {
            panic!("❌ 输入文件格式错误：列数不足。请检查是否是经过 clean_vcf 处理的文件。");
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
    /// 从 VCF 文件的一行解析出 SnpInfo
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

        // 解析 format 和 sample 数据
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

    /// 输出为一行 TSV 格式
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
    /// 创建对象
    pub fn from_path(path: &str) -> Self {
        Self {
            summary_path: path.to_string(),
        }
    }
    /// 解析 VCF 文件，返回清洗后的 SNP 列表
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

    /// 从 cleaned GWAS summary 文件加载所有 SNP 信息
    pub fn load_snps(&self) -> io::Result<Vec<SnpInfo>> {
        let file = File::open(&self.summary_path)?;
        let reader = BufReader::new(file);
        let mut snps = Vec::new();
        for (_i, line) in reader.lines().enumerate() {
            let line = line?;
            // 跳过空行或表头
            if line.trim().is_empty() {
                continue;
            }
            if line.contains("CHR") {
                continue; // 跳过表头行
            }

            if let Some(snp_info) = SnpInfo::from_cleaned_line(&line) {
                snps.push(snp_info);
            }
        }

        Ok(snps)
    }

    /// 获取文件名（不带路径）
    pub fn filename(&self) -> Option<String> {
        std::path::Path::new(&self.summary_path)
            .file_name()
            .and_then(|s| s.to_str())
            .map(|s| s.to_string())
    }

    /// 过滤 p 值小于阈值
    pub fn filter_by_p(snps: Vec<SnpInfo>, p_threshold: f64) -> Vec<SnpInfo> {
        let result: Vec<SnpInfo> = snps
            .into_iter()
            .filter(|s| s.p.map_or(false, |p| p <= p_threshold))
            .collect();
        println!(
            "filter by p value<{},remaining snp {}",
            p_threshold,
            result.len()
        );
        result
    }

    /// 过滤 p 值小于阈值
    pub fn filter_by_maf(snps: Vec<SnpInfo>, maf_threshold: f64) -> Vec<SnpInfo> {
        let result: Vec<SnpInfo> = snps
            .into_iter()
            .filter(|s| s.maf.map_or(false, |maf| maf >= maf_threshold))
            .collect();
        println!(
            "filter by maf>{},remaining snp {}",
            maf_threshold,
            result.len()
        );
        result
    }

    /// 过滤 SNP 名称列表
    pub fn merge_alleles<'a>(
        snps: Vec<SnpInfo>,
        snp_list: Vec<(&'a str, &'a str, &'a str)>,
    ) -> Vec<SnpInfo> {
        let set: HashSet<(&str, &str, &str)> = snp_list.into_iter().collect();
        let mut result1: Vec<SnpInfo> = snps
            .clone()
            .into_iter()
            .filter(|s| {
                let key = (s.snp.as_str(), s.a2.as_str(), s.a1.as_str());
                set.contains(&key)
            })
            .map(|mut s| {
                if let Some(maf) = s.maf {
                    s.maf = Some(1.0 - maf);
                }
                s
            })
            .collect();
        let result2: Vec<SnpInfo> = snps
            .into_iter()
            .filter(|s| {
                let key = (s.snp.as_str(), s.a1.as_str(), s.a2.as_str());
                set.contains(&key)
            })
            .collect();
        result1.extend(result2);
        println!("merge alleles, remaining snps {}", result1.len());
        result1
    }

    pub fn compute_h2(&self, ld: LdFile) -> io::Result<f64> {
        // 1️⃣ 加载 GWAS summary 中的 SNP
        let snps = self.load_snps()?; // Vec<SnpInfo>

        // 2️⃣ 加载 LD 分数
        let ld_map = ld.load_ld_scores()?; // HashMap<snp_name, ldscore>

        // 3️⃣ 遍历 SNP，计算每个 SNP 对 h² 的贡献
        let mut x_vec = Vec::new(); // L2
        let mut y_vec = Vec::new(); // Z^2 -1
        let mut neff_vec = Vec::new(); // NEFF

        for snp in snps.iter() {
            // 取 beta
            let beta = match snp.beta {
                Some(b) => b,
                None => continue, // 没有 beta 的 SNP 跳过
            };
            let se = match snp.se {
                Some(b) => b,
                None => continue, // 没有 beta 的 SNP 跳过
            };
            let neff = match snp.neff {
                Some(b) => b,
                None => continue,
            };
            let z = beta / se;

            let l2 = match ld_map.get(&snp.snp) {
                Some(&val) => val,
                None => continue, // LD 不存在则跳过
            };

            x_vec.push(l2);
            y_vec.push(z * z - 1.0);
            neff_vec.push(neff);
        }
        let m = x_vec.len() as f64; // SNP 总数
        let mean_neff = neff_vec.iter().sum::<f64>() / m;

        // 最小二乘法计算 slope
        let x_mean = x_vec.iter().sum::<f64>() / m;
        let y_mean = y_vec.iter().sum::<f64>() / m;

        let numerator: f64 = x_vec
            .iter()
            .zip(y_vec.iter())
            .map(|(x, y)| (x - x_mean) * (y - y_mean))
            .sum();

        let denominator: f64 = x_vec.iter().map(|x| (x - x_mean).powi(2)).sum();
        let slope = numerator / denominator;

        // h^2 = slope * M / N
        let h2 = slope * m / mean_neff;
        println!(
            "slope {},snp number {},mean effective sample number {}",
            slope, m, mean_neff
        );
        Ok(h2)
    }
    /// 将 SNP 列表写入任意 Writer
    pub fn write_snps(snps: &[SnpInfo], output_path: &str) -> io::Result<()> {
        let output_file = File::create(output_path)?;
        let mut writer = BufWriter::new(output_file);
        writeln!(writer, "CHR\tPOS\tSNP\tA1\tA2\tBETA\tSE\tP\tMAF\tNEFF")?;
        for snp in snps {
            writeln!(writer, "{}", snp.to_tsv())?;
        }
        Ok(())
    }
}

impl LdFile {
    /// 创建对象
    pub fn from_path(path: &str) -> Self {
        Self {
            ld_path: path.to_string(),
        }
    }
    /// 从 LD Score 文件中读取 SNP -> LD Score 的映射
    pub fn load_ld_scores(&self) -> io::Result<HashMap<String, f64>> {
        let file = File::open(&self.ld_path)?;
        let reader = BufReader::new(file);
        let mut ld_map = HashMap::new();
        // 跳过表头
        for (_i, line) in reader.lines().enumerate() {
            let line = line?;
            if line.contains("SNP") {
                continue;
            }
            let fields: Vec<&str> = line.split_whitespace().collect();
            if fields.len() < 6 {
                continue; // 跳过不完整行
            }
            let snp = fields[1].to_string(); // SNP 列
            let ld_score = fields[5].parse::<f64>().unwrap_or(0.0); // L2 列
            ld_map.insert(snp, ld_score);
        }

        Ok(ld_map)
    }
}

impl RefAlleles {
    pub fn from_path(path: &str) -> Self {
        Self {
            alleles_path: path.to_string(),
        }
    }
    pub fn load_alleles(&self) -> io::Result<Vec<AllelesInfo>> {
        let file = File::open(&self.alleles_path)?;
        let reader = BufReader::new(file);
        let mut snps = Vec::new();
        for (_i, line) in reader.lines().enumerate() {
            let line = line?;
            // 跳过空行或表头
            if line.trim().is_empty() {
                continue;
            }
            if line.contains("CHR") {
                continue; // 跳过表头行
            }
            if let Some(snp_info) = AllelesInfo::from_cleaned_line(&line) {
                snps.push(snp_info);
            }
        }
        Ok(snps)
    }
    pub fn extract_snp_names(snps: &[AllelesInfo]) -> Vec<(&str, &str, &str)> {
        snps.iter()
            .map(|s| (s.snp.as_str(), s.a2.as_str(), s.a1.as_str()))
            .collect()
    }
}

impl AllelesInfo {
    pub fn from_cleaned_line(line: &str) -> Option<Self> {
        if line.starts_with('#') || line.trim().is_empty() {
            return None;
        }
        let fields: Vec<&str> = line.trim().split('\t').collect();

        if fields.len() < 3 {
            panic!("❌ 输入文件格式错误：列数不足。请检查是否是经过 clean_vcf 处理的文件。");
        }

        Some(Self {
            snp: fields[0].to_string(),
            a1: fields[1].to_string(),
            a2: fields[2].to_string(),
        })
    }
}

/// compute the genetic correlation
pub fn compute_genetic_correlation(
    file_path1: GwasSummary,
    file_path2: GwasSummary,
    ld_path: LdFile,
) -> io::Result<f64> {
    // 计算两个性状的 SNP 遗传力
    let h1_2 = file_path1.compute_h2(ld_path.clone())?;
    let h2_2 = file_path2.compute_h2(ld_path.clone())?;
    let ld_map = ld_path.load_ld_scores()?;
    let h1_h2 = (h1_2 * h2_2).sqrt();
    let mut map1_nol: HashMap<String, (f64, f64)> = HashMap::new();
    let mut map2_nol: HashMap<String, (f64, f64)> = HashMap::new();

    let snps1 = file_path1.load_snps()?;
    let snps2 = file_path2.load_snps()?;
    for snp in snps1.iter() {
        // 取 beta
        let beta = match snp.beta {
            Some(b) => b,
            None => continue, // 没有 beta 的 SNP 跳过
        };
        let se = match snp.se {
            Some(b) => b,
            None => continue, // 没有 beta 的 SNP 跳过
        };
        let neff = match snp.neff {
            Some(b) => b,
            None => continue,
        };
        let z = beta / se;

        map1_nol.insert(snp.snp.clone(), (z, neff));
    }

    for snp in snps2.iter() {
        // 取 beta
        let beta = match snp.beta {
            Some(b) => b,
            None => continue, // 没有 beta 的 SNP 跳过
        };
        let se = match snp.se {
            Some(b) => b,
            None => continue, // 没有 beta 的 SNP 跳过
        };
        let neff = match snp.neff {
            Some(b) => b,
            None => continue,
        };
        let z = beta / se;

        map2_nol.insert(snp.snp.clone(), (z, neff));
    }

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

    // 构建 Z1Z2 和 L2 向量
    let mut z1z2_vec: Vec<f64> = Vec::new();
    let mut l2_vec: Vec<f64> = Vec::new();

    for key in &inter_keys {
        let (z1, l2_1, _) = map1_filtered.get(key).unwrap();
        let (z2, l2_2, _) = map2_filtered.get(key).unwrap();
        z1z2_vec.push(z1 * z2);
        // 两个文件的 L2 理论上应该一样，这里取平均更稳妥
        l2_vec.push((l2_1 + l2_2) / 2.0);
    }

    // 最小二乘法线性回归 slope
    let x_mean = l2_vec.iter().sum::<f64>() / m;
    let y_mean = z1z2_vec.iter().sum::<f64>() / m;

    let numerator: f64 = l2_vec
        .iter()
        .zip(z1z2_vec.iter())
        .map(|(x, y)| (x - x_mean) * (y - y_mean))
        .sum();

    let denominator: f64 = l2_vec.iter().map(|x| (x - x_mean).powi(2)).sum();

    let slope = numerator / denominator;

    // 根据理论公式计算 genetic correlation
    let genetic_correlation = slope * m / ((n1 * n2).sqrt() * h1_h2);

    println!(
        "M = {:.0}, N1 = {:.0}, N2 = {:.0}, h1² = {:.4}, h2² = {:.4}, slope = {:.6}, r_g = {:.6}",
        m, n1, n2, h1_2, h2_2, slope, genetic_correlation
    );

    Ok(genetic_correlation)
}
