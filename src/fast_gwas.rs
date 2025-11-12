use std::{
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    path::Path,
    time::Instant,
};
use rustc_hash::FxHashMap;

/// 筛选 GWAS 文件，只保留 HapMap3 中的 SNP，并处理等位基因方向与 MAF
pub fn filter_gwas(hm3_path: &Path, gwas_path: &Path, out_path: &Path) -> std::io::Result<()> {
    let start = Instant::now();

    // === 1️⃣ 加载 HapMap3 SNPs ===
    let mut hm3_map: FxHashMap<String, (String, String)> = FxHashMap::default();
    let hm3_file = BufReader::with_capacity(16 * 1024 * 1024, File::open(hm3_path)?);
    for (i, line) in hm3_file.lines().enumerate() {
        let l = line?;
        if i == 0 || l.is_empty() {
            continue;
        }
        let mut iter = l.split_whitespace();
        if let (Some(snp), Some(a1), Some(a2)) = (iter.next(), iter.next(), iter.next()) {
            hm3_map.insert(snp.to_owned(), (a1.to_owned(), a2.to_owned()));
        }
    }
    println!("Loaded {} HapMap3 SNPs", hm3_map.len());

    // === 2️⃣ 筛选 GWAS 文件 ===
    let gwas_file = BufReader::with_capacity(64 * 1024 * 1024, File::open(gwas_path)?);
    let mut writer = BufWriter::with_capacity(32 * 1024 * 1024, File::create(out_path)?);
    let mut total = 0usize;
    let mut kept = 0usize;
    let mut header_written = false;

    for line in gwas_file.lines() {
        let line = line?;
        if !header_written {
            writeln!(writer, "{}", line)?;
            header_written = true;
            continue;
        }
        total += 1;
        if line.is_empty() {
            continue;
        }

        // 手动分割以减少 Vec 分配
        let mut parts = line.split('\t');
        let (chr, pos, snp, a1, a2, beta_s, se, p, maf_s, neff) = match (
            parts.next(),
            parts.next(),
            parts.next(),
            parts.next(),
            parts.next(),
            parts.next(),
            parts.next(),
            parts.next(),
            parts.next(),
            parts.next(),
        ) {
            (Some(c), Some(p1), Some(s), Some(a_1), Some(a_2), Some(b), Some(se_), Some(p_), Some(m), Some(n)) => {
                (c, p1, s, a_1, a_2, b, se_, p_, m, n)
            }
            _ => continue,
        };

        let Some((hm3_a1, hm3_a2)) = hm3_map.get(snp) else { continue };

        let Ok(mut beta) = beta_s.parse::<f64>() else { continue };
        let Ok(mut maf) = maf_s.parse::<f64>() else { continue };

        if a1 == hm3_a1 && a2 == hm3_a2 {
            // ok
        } else if a1 == hm3_a2 && a2 == hm3_a1 {
            beta = -beta;
            maf = 1.0 - maf;
        } else {
            continue;
        }

        if maf <= 0.01 {
            continue;
        }

        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{:.6}\t{}\t{}\t{:.6}\t{}",
            chr, pos, snp, hm3_a1, hm3_a2, beta, se, p, maf, neff
        )?;
        kept += 1;
    }

    writer.flush()?;
    println!("✅ Done! Filtered {} / {} lines", kept, total);
    println!("⏱️ 程序运行时间: {:.2?}", start.elapsed());
    Ok(())
}