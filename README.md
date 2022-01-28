# Strawberry_chiloensis_GP33

This repository contains the scripts used to produce the figures of Cauret et al. 2022: Chromosome-scale assembly with a phased sex-determining region resolves features of early Z and W chromosome differentiation in a wild octoploid strawberry

`Fst_last.R`: produces the supplemental Fst figures

`cov.R`: produces the Fig. 3, i.e. combination of a dotplot of the Z and W sequence, coverage analysis (difference between females and males), and repeat coverage.

`all_sex_chr_analysis_last.R` produces Fig. 5, i.e. repeat coverage, intersexual weighted FST and difference between males and females in coverage (log2(mean F) - log2(mean M)) across the sex chromosome

`dotplot.R`: produces 4 dotplots (3 genome wide, 1 subset of homeologs of group 6): figures S1, S2, S3 and main Fig. 2 

Input files:

- FST files were produced with VCFtools (0.1.17) 

```python3 /dfs/Liston_Lab/workspace/caroline/scripts/run_vcftools_pop_analysis.py --vcf_files /dfs/Liston_Lab/scratch/cauretc/2021_sex_chrom_analysis/wgs/vcf_files/*_filtered_ind_filtered_qual.vcf.gz --pop1_file ../13_F_M/13_F.txt --pop2_file ../13_F_M/13_M.txt --out_dir /dfs/Liston_Lab/scratch/cauretc/2021_sex_chrom_analysis/wgs/pop_stats/13_F_M_no_miss_filt```

`run_vcftools_pop_analysis.py` simply runs the command `vcftools --gzvcf [file] --weir-fst-pop [pop1] --weir-fst-pop [pop2] --out [out_prefix] --fst-window-size [size]` in an autommated way in the different files conatains in a directory. 

- coverage files were produced with samtools 1.14 (Using htslib 1.14)

```samtools depth -H -a -f /dfs/Liston_Lab/scratch/cauretc/2021_sex_chrom_analysis/wgs/bam_files/Liston_bam_list_chil.txt -o /dfs/Liston_Lab/scratch/cauretc/2021_sex_chrom_analysis/wgs/depth_files/no_secondary/Liston_Fvb6-1_depth_no_qual_no_sec.txt -r Fvb6-1_RagTag```
