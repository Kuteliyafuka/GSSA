# troduction


This project aims to develop a flexible tool to perform genetic analyses using GWAS summary statistics.  

The main reasons for creating our own tool are:

1. The analysis requires multiple steps and tools, but the summary statistics themselves are not very complex.
2. The corresponding analyses (heritability, genetic correlation, MR, colocalization, etc.) can be implemented efficiently with a unified framework.

Below is a summary of the main tasks during my second rotation:

## Task Overview

| Mission                                      | Tool         | Method                     |
| ------------------------------------------- | ------------ | -------------------------- |
| Compute heritability and genetic correlation | LDSC         | Linear regression          |
| Choose instrument SNPs                       | PLINK        |                            |
| Infer causal inference                       | TwoSampleMR  | Weighted linear regression |
| Identify pleiotropic loci                    | ASSET / pleio |                            |
| Identify causal variant                      | coloc        | Bayesian                   |

---

## Workflow and Implementation

This workflow illustrates the main steps using two example datasets:

- `ieu-a-2`: BMI dataset (exposure)
- `ieu-a-7`: CHD dataset (outcome)

### 1. Compute Heritability and Genetic Correlation

```bash
# Step 1: Clean the VCF file
gwas_tool --vcf-clean ieu-a-2.vcf --output ieu-a-2.sumstats

# Step 2: Flip alleles and filter by MAF
gwas_tool --modify ieu-a-2.sumstats \
          --merge-alleles w_hm3.snplist.txt \
          --min-maf 0.01 \
          --output ieu-a-2_modified.sumstats

# Step 3: Compute SNP heritability
gwas_tool --compute-h2 ieu-a-2_modified.sumstats --ld summary.ldscore

# Step 4: Compute genetic correlation
gwas_tool --genetic-correlation ieu-a-2_modified.sumstats \
          ieu-a-7_modified.sumstats
````

> **Note:** `summary.ldscore` is generated from the European 1000 Genomes panel: `cat eur_w_ld_chr/*.gz > summary.ldscore.gz`.

---

### 2. Choose Instrument SNPs and Perform MR

```bash
# Keep SNPs with P < 5e-8 in exposure
gwas_tool --modify ieu-a-2_modified.sumstats \
          --max-p 5e-8 \
          --output ieu-a-2_modified2.sumstats

# SNP clumping (recommended to use PLINK for speed)
gwas_tool --snp-clump ieu-a-2_modified2.sumstats \
          --ld-dir eur_ld \
          --output ieu-a-2_clumped.sumstats

# Filter outcome by clumped SNPs
gwas_tool --modify ieu-a-7_modified.sumstats \
          --filter-by-snps ieu-a-2_clumped.sumstats \
          --output ieu-a-7_clumped.sumstats

# Perform MR analysis
gwas_tool --mr-ivm --exposure ieu-a-2_clumped.sumstats \
          --outcome ieu-a-7_clumped.sumstats

gwas_tool --mr-egger --exposure ieu-a-2_clumped.sumstats \
           --outcome ieu-a-7_clumped.sumstats
```

---

### 3. Identify Pleiotropic Loci

```bash
gssa --identify pleio-identify BMI_filtered.sum MMD_filtered.sum \
     --min-z 8 --max-p 5e-8 \
     --output pleiotropic_loci.tsv

# Filter MMD dataset by pleiotropic SNPs
gssa --modify MMD_filtered.sum \
     --filter-by-snps pleiotropic_loci.tsv \
     --output MMD_pleio.sum
```

> Optional: Run PLINK SNP clumping on `MMD_pleio.sum` to generate clumped SNPs for downstream analysis.

---

### 4. Identify Lead SNPs and Perform Colocalization

```bash
# Filter MMD dataset by clumped SNPs
gssa --modify MMD_filtered.sum \
     --filter-by-snps snps.clumped \
     --output lead_snps.txt

# Colocalization between MDD and BMI
gssa --lead-snps lead_snps.txt \
     --coloc BMI_filtered.sum MMD_filtered.sum \
     --output coloc_result.tsv
```

---

### 5. Fast Filtering Using Default Parameters

For improved performance, `fast-modify` can quickly filter summary statistics using default parameters:

```bash
gssa --fast-modify ieu-a-2.sumstats \
     --merge-alleles w_hm3.snplist.txt
```

---

## Summary

This workflow demonstrates the core functionality of the tool:

1. Data cleaning and SNP filtering.
2. Heritability and genetic correlation estimation.
3. Instrument SNP selection and MR analysis.
4. Pleiotropic locus identification.
5. Colocalization analysis.

With these modules, users can perform end-to-end GWAS summary statistic analyses using a unified framework.

