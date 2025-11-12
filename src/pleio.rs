use std::fs::File;
use std::io::{BufWriter, Write};
use crate::gwas_sum::SnpInfo;
use std::collections::HashMap;
#[derive(Debug, Clone)]
pub struct PleioResult {
    pub snp: String,
    pub z1: f64,
    pub z2: f64,
    pub p1: f64,
    pub p2: f64,
    pub abs_sum: f64,
    pub abs_joint: f64,
    pub pleio_type: String, // "concordant", "opposite", "non"
}
#[derive(Debug)]
pub struct ColocResult {
    pub lead_snp: String,
    pub pleio_region: String,
    pub pp_h0: f64,
    pub pp_h1: f64,
    pub pp_h2: f64,
    pub pp_h3: f64,
    pub pp_h4: f64,
    pub max_h4_snp: String,
    pub max_h4_value: f64,
}

/// 主函数：识别多效性 SNP
pub fn pleio_identify(snps1: Vec<SnpInfo>, snps2: Vec<SnpInfo>, z_threshold: f64,p_threshold:f64) -> Vec<PleioResult> {
    assert_eq!(
        snps1.len(),
        snps2.len(),
        "Input vectors must have the same length"
    );

    let mut results = Vec::new();

    for (s1, s2) in snps1.iter().zip(snps2.iter()) {
        let (beta1, se1, p1) = match (s1.beta, s1.se, s1.p) {
            (Some(b), Some(s), Some(p)) => (b, s, p),
            _ => continue,
        };
        let (beta2, se2, p2) = match (s2.beta, s2.se, s2.p) {
            (Some(b), Some(s), Some(p)) => (b, s, p),
            _ => continue,
        };

        // 过滤条件
        if p1 >= p_threshold || p2 >= p_threshold {
            continue;
        }

        let z1 = beta1 / se1;
        let z2 = beta2 / se2;
        let abs_sum = z1.abs() + z2.abs();
        let abs_joint = (z1 + z2).abs();

        let pleio_type = if abs_sum > z_threshold || abs_joint > z_threshold {
            if (abs_sum - abs_joint).abs() < 1e-6 && abs_joint > z_threshold {
                "concordant"
            } else if abs_sum > abs_joint {
                "opposite"
            } else {
                "non"
            }
        } else {
            "non"
        }
        .to_string();

        results.push(PleioResult {
            snp: s1.snp.clone(),
            z1,
            z2,
            p1,
            p2,
            abs_sum,
            abs_joint,
            pleio_type,
        });
    }

    results
}

/// 写出结果到三个文件
pub fn write_pleio_results(results: &[PleioResult], output_path: &str) -> std::io::Result<()> {
    let base = if let Some(pos) = output_path.rfind('.') {
        &output_path[..pos]
    } else {
        output_path
    };

    // 构造三个输出文件名
    let cons_path = format!("{}_cons_pleio.txt", base);
    let oppo_path = format!("{}_oppo_pleio.txt", base);
    let non_path = format!("{}_non_pleio.txt", base);

    // 打开写入器
    let mut cons = BufWriter::new(File::create(&cons_path)?);
    let mut oppo = BufWriter::new(File::create(&oppo_path)?);
    let mut non = BufWriter::new(File::create(&non_path)?);

    let header = "SNP\tz1\tz2\tp1\tp2\t|z1+z2|\t|z1|+|z2|\n";
    for f in [&mut cons, &mut oppo, &mut non] {
        f.write_all(header.as_bytes())?;
    }

    for r in results {
        let line = format!(
            "{}\t{:.4}\t{:.4}\t{:.3e}\t{:.3e}\t{:.4}\t{:.4}\n",
            r.snp, r.z1, r.z2, r.p1, r.p2, r.abs_joint, r.abs_sum
        );
        match r.pleio_type.as_str() {
            "concordant" => cons.write_all(line.as_bytes())?,
            "opposite" => oppo.write_all(line.as_bytes())?,
            _ => non.write_all(line.as_bytes())?,
        }
    }

    Ok(())
}

pub fn coloc(snps1: Vec<SnpInfo>, snps2: Vec<SnpInfo>, lead_snp: &str, output_path: String) -> std::io::Result<ColocResult> {

    // --- 参数设置 ---
    let w = 0.04_f64;   
    let p1 = 1e-4_f64;
    let p2 = 1e-4_f64;
    let p12 = 1e-5_f64;

    // --- 建立SNP索引 ---
    let mut map2: HashMap<String, &SnpInfo> = HashMap::new();
    for s in &snps2 {
        map2.insert(s.snp.clone(), s);
    }

    // --- 计算ABF ---
    let abf_fun = |beta: f64, v: f64, w: f64| -> f64 {
        (v / (v + w)).sqrt() * ((w * beta.powi(2)) / (2.0 * v * (v + w))).exp()
    };

    let mut abf1 = Vec::new();
    let mut abf2 = Vec::new();
    let mut snps_common = Vec::new();

    for s1 in &snps1 {
        if let Some(s2) = map2.get(&s1.snp) {
            if let (Some(b1), Some(se1), Some(b2), Some(se2)) = (s1.beta, s1.se, s2.beta, s2.se) {
                let v1 = se1.powi(2);
                let v2 = se2.powi(2);
                abf1.push(abf_fun(b1, v1, w));
                abf2.push(abf_fun(b2, v2, w));
                snps_common.push(s1.snp.clone());
            }
        }
    }

    if abf1.is_empty() {
        eprintln!("⚠️ No overlapping SNPs for lead SNP: {}", lead_snp);
        return Ok(ColocResult {
            lead_snp: lead_snp.to_string(),
            pleio_region: "NA".to_string(),
            pp_h0: f64::NAN,
            pp_h1: f64::NAN,
            pp_h2: f64::NAN,
            pp_h3: f64::NAN,
            pp_h4: f64::NAN,
            max_h4_snp: "NA".to_string(),
            max_h4_value: f64::NAN,
        });
    }

    // --- 汇总ABF ---
    let sum_abf1: f64 = abf1.iter().sum();
    let sum_abf2: f64 = abf2.iter().sum();
    let sum_prod: f64 = abf1.iter().zip(&abf2).map(|(x, y)| x * y).sum();

    // --- 五个假设的贝叶斯因子 ---
    let bf_h0 = 1.0;
    let bf_h1 = p1 * sum_abf1;
    let bf_h2 = p2 * sum_abf2;
    let bf_h3 = p1 * p2 * (sum_abf1 * sum_abf2 - sum_prod);
    let bf_h4 = p12 * sum_prod;

    let bfs = vec![bf_h0, bf_h1, bf_h2, bf_h3, bf_h4];
    let sum_bf: f64 = bfs.iter().sum();

    // --- 后验概率 ---
    let pp: Vec<f64> = bfs.iter().map(|x| x / sum_bf).collect();
    let pp_h0 = pp[0];
    let pp_h1 = pp[1];
    let pp_h2 = pp[2];
    let pp_h3 = pp[3];
    let pp_h4 = pp[4];

    // --- 每个SNP的H4概率 ---
    let snp_pp_h4: Vec<f64> = abf1.iter().zip(&abf2).map(|(x, y)| x * y / sum_prod).collect();

    // 找出 H4 最高的 SNP
    let (max_idx, max_h4_value) = snp_pp_h4
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
        .unwrap();
    let max_h4_snp = snps_common[max_idx].clone();
    let pleio_region = snps1.iter().find(|s| s.snp == lead_snp).map_or("NA".to_string(), |s| format!("chr{}:{}-{}", s.chr, s.pos.saturating_sub(500_000), s.pos + 500_000));
    // 写出单 locus 文件
    let file = File::create(&output_path)?;
    let mut writer = BufWriter::new(file);
    writeln!(writer, "SNP\tPP_H4")?;
    for (snp, pp) in snps_common.iter().zip(snp_pp_h4.iter()) {
        writeln!(writer, "{}\t{:.6e}", snp, pp)?;
    }
    writer.flush()?;

    Ok(ColocResult {
        lead_snp: lead_snp.to_string(),
        pleio_region,
        pp_h0,
        pp_h1,
        pp_h2,
        pp_h3,
        pp_h4,
        max_h4_snp,
        max_h4_value: *max_h4_value,
    })
}



pub fn multi_coloc(
    lead_snps: Vec<SnpInfo>,
    snps1: Vec<SnpInfo>,
    snps2: Vec<SnpInfo>,
    output_path: String,
) -> std::io::Result<()> {

    let mut results = Vec::new();

    for lead in lead_snps {
        let snp1_subset: Vec<SnpInfo> = snps1
            .iter()
            .filter(|s| s.chr == lead.chr && (s.pos as i64 - lead.pos as i64).abs() <= 250_000)
            .cloned()
            .collect();
        let snp2_subset: Vec<SnpInfo> = snps2
            .iter()
            .filter(|s| s.chr == lead.chr && (s.pos as i64 - lead.pos as i64).abs() <= 250_000)
            .cloned()
            .collect();

        if snp1_subset.is_empty() || snp2_subset.is_empty() {
            eprintln!("⚠️ Skipping lead SNP {} (no overlapping region)", lead.snp);
            continue;
        }

        let coloc_output_path =
            format!("{}_coloc_{}.txt", output_path.trim_end_matches(".txt"), lead.snp);

        match coloc(snp1_subset, snp2_subset, &lead.snp, coloc_output_path) {
            Ok(res) => results.push(res),
            Err(err) => eprintln!("Error in coloc for {}: {:?}", lead.snp, err),
        }
    }

    // --- 输出汇总结果 ---
    let file = File::create(&output_path)?;
    let mut writer = BufWriter::new(file);
    writeln!(writer, "SNP\tpleio_region\tPP_H0\tPP_H1\tPP_H2\tPP_H3\tPP_H4\tmax_H4_snp\tmax_H4_value\tcausal_variant")?;
    let mut causal_variants_number = 0;
    let mut is_causal_variant: &str;
    for r in results {
        if r.pp_h4 > 0.8 {
            is_causal_variant = "True";
            causal_variants_number += 1;
        } else {
            is_causal_variant = "False";
        }
        writeln!(
            writer,
            "{}\t{}\t{:.6e}\t{:.6e}\t{:.6e}\t{:.6e}\t{:.6e}\t{}\t{:.6e}\t{}",
            r.lead_snp, r.pleio_region,r.pp_h0,r.pp_h1, r.pp_h2, r.pp_h3, r.pp_h4, r.max_h4_snp, r.max_h4_value, is_causal_variant
        )?;
    }
    writer.flush()?;
    println!("✅ Coloc analysis completed. causal variant number {}", causal_variants_number);
    println!("✅ Summary results written to: {}", output_path);
    Ok(())
}
