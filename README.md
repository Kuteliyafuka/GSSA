# 阶段性小结



## Motivation

为什么要自己写一个工具？

1、因为这一部分的分析需要使用若干工具，而且gwas summary statistics数据本身不复杂，所以相应的分析也并不复杂，完全可以自己实现。

一下是我二轮轮转期间的主要任务

## 任务清单

| mission                                      | tool         | method                     |
| -------------------------------------------- | ------------ | -------------------------- |
| compute heritability and genetic correlation | LDSC         | linear regression          |
| choose instrument SNPs                       | Plink        |                            |
| infer causal inference                       | TwoSamplesMR | weighted linear regression |
| identify pleiotropic loci                    | asset/pleio  |                            |
| identify causal vaciance                     | coloc        | bayes                      |



## 具体实现

首先，下载ieu-a-2数据集和ieu-a-7数据集。

ieu-a-2是BMI的数据集，ieu-a-7是CHD（冠心病）数据集

```bash
# Part1: compute h2 and genetic correlation
# clean the vcf data into tsv format
gwas_tool --vcf-clean ieu-a-2.vcf --output ieu-a-2.sumstats

# flip alleles and filter MAF<0.01
gwas_tool --modify ieu-a-2.sumstats --merge_alleles w_hm3.snplist.txt --min-maf 0.01 --output ieu-a-2_modified.sumstats

# compute h2
# note: summary.ldscore is got by in eur_w_ld_chr(the european 1000 genome panel), cat *.gz > summary.ldscore.gz
gwas_tool --compute-h2 ieu-a-2_modified.sumstats --ld summary.ldscore

# compute genetic correlation
gwas_tool --genetic-correlation ieu-a-2_modified.sumstats ieu-a-7_modified.sumstats



# Part2: choose instrument SNPs and do MR
# remain the snps with P < 5e-8 in exposure data
gwas_tool --modify ieu-a-2_modified.sumstats --max-p 5e-8 --output ieu-a-2_modified2.sumstats
# snp clump step is strongly recommended to run with plink, even gwas_tool has snp_clump function, it is too slow.
# first, use plink to generate eur_ld.ld file(which stores ld r2 data),then split eur_ld.ld into eur_ld/chr?.ld files.
gwas_tool --snp-clump ieu-a-2_modified2.sumstats --ld-dir eur_ld --output ieu-a-2_clumped.sumstats

# filter the outcome data with the clumped snps
gwas_tool --modify ieu-a-7_modified.sumstats --filter-by-snps ieu-a-2_clumped.sumstats --output ieu-a-7_clumped.sumstats

# do the mr-ivm and mr-egger
gwas_tool --mr-ivm --exposure ieu-a-2_clumped.sumstats --outcome ieu-a-7_clumped.sumstats
gwas_tool --mr-egger --exposure ieu-a-2_clumped.sumstats --outcome ieu-a-7_clumped.sumstats
```

问题1：这里的h2的计算和ldsc的计算有点差别，可能是没有设置权重的原因（用同一个文件），这个问题还没有解决

问题2：--snp-clump这一步速度太慢了，依然还是需要用plink来做snp clump
