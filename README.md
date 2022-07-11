# Strawberry_chiloensis_GP33

This repository contains the scripts used to produce the figures of Cauret et al. 2022: Chromosome-scale assembly with a phased sex-determining region resolves features of early Z and W chromosome differentiation in a wild octoploid strawberry

- `Fst_last.R`: produces the supplemental Fst figures

- `cov.R`: produces the Fig. 3, i.e. combination of a dotplot of the Z and W sequence, coverage analysis (difference between females and males using normalized coverage files), and repeat coverage.

- `all_sex_chr_analysis_last.R` produces Fig. 5, i.e. repeat coverage, intersexual weighted FST and difference between males and females in read coverage (log2(mean F) - log2(mean M)) across the sex chromosome

- `dotplot.R`: produces 4 dotplots (3 genome wide, 1 subset of homeologs of group 6): figures S1, S2, S3 and main Fig. 2 

- `phylo_SDR_gene.R`: produces a phylogeny figure combining three gene trees: main figure 6. 

Input files:

- FST files were produced with VCFtools (0.1.17) 

```python3 /dfs/Liston_Lab/workspace/caroline/scripts/run_vcftools_pop_analysis.py --vcf_files /dfs/Liston_Lab/scratch/cauretc/2021_sex_chrom_analysis/wgs/vcf_files/*_filtered_ind_filtered_qual.vcf.gz --pop1_file ../13_F_M/13_F.txt --pop2_file ../13_F_M/13_M.txt --out_dir /dfs/Liston_Lab/scratch/cauretc/2021_sex_chrom_analysis/wgs/pop_stats/13_F_M_no_miss_filt```

`run_vcftools_pop_analysis.py` simply runs the command `vcftools --gzvcf [file] --weir-fst-pop [pop1] --weir-fst-pop [pop2] --out [out_prefix] --fst-window-size [size]` in an autommated way in the different files conatains in a directory. 

- raw coverage files were produced with samtools 1.14 (Using htslib 1.14)

```samtools depth -H -a -f /dfs/Liston_Lab/scratch/cauretc/2021_sex_chrom_analysis/wgs/bam_files/Liston_bam_list_chil.txt -o /dfs/Liston_Lab/scratch/cauretc/2021_sex_chrom_analysis/wgs/depth_files/no_secondary/Liston_Fvb6-1_depth_no_qual_no_sec.txt -r Fvb6-1_RagTag```

- `normalization.py` (courtesy of Dr. Andrew S. Tupper) was used on the raw coverage files to obtain the normalized coverage (data normalized by the coverage of a representative autosome) per individual, obtain the mean ratio (female/male) of the coverage and calculate non-overlapping window average. The final output file (`ratio_merged_window.txt`) was used as input for the above R scripts.
```
# Calculate the means from the autosomal files
python3 normalization.py means Liston_Fvb3-1_depth_no_qual_no_sec.txt UCD_Fvb3-1_depth_no_qual_no_sec.txt #produces means.db
# Use previously calculated means to normalize the data, then use the female and male ID files to calculate the mean coverage per sex and determine the female/male ratio at each position
# The ID files contain 1 individual per line
python3 normalization.py ratio data/Males_ID_cov.txt data/Females_ID_cov.txt data/*3-1*.txt > ratio_3_1.txt #uses means.db without having to specify it
python3 normalization.py ratio data/Males_ID_cov.txt data/Females_ID_cov.txt data/*6-1*.txt > ratio_6_1.txt
python3 normalization.py ratio data/Males_ID_cov.txt data/Females_ID_cov.txt data/*atg*.txt > ratio_atg.txt
# Calculate averages from ratio files on 10kb non-overlapping windows
python3 normalization.py window 10000 ratio_3_1.txt ratio_6_1.txt ratio_atg.txt > ratio_merged_window.txt
```

## GDR submission

The final genome assembly is available on GDR (https://www.rosaceae.org/).

Since the chromosome nomenclature was decided at the end of the analysis, the naming for the figure was performed within the scripts described above. For GDR, I made a script to rename the genome assembly and annotation with the final nomenclature (i.e. "Fchil1-B1" instead of "Fchil1-B1"). 
Run: 
```
python3 rename_genome.py --nomenclature_file chromosome_nomenclature.txt --fasta ~/Documents/Caroline_postdoc_OSU/2022_G3_zenodo_temp/ragtag.scaffolds.curated.reorientated.fasta --gff_file ~/Documents/Caroline_postdoc_OSU/2022_G3_zenodo_temp/Fchil_GP33_withoutW_v0.1.gff_polished --output_fasta Fchil_GP33_v1.0.fa --output_gff_file Fchil_GP33_v1.0.gff
```