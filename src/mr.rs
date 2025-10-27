use crate::gwas_sum::SnpInfo;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};

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
        if i == 0 { continue; } // 跳过表头

        let cols: Vec<&str> = line.split_whitespace().collect();
        if cols.len() < 7 { continue; }

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
                println!("{}",col_name);
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