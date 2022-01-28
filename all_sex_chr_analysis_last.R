library("dplyr")
library("plyr")
library("ggplot2")
library("stringr")
library("grid")
require("data.table")
library("cowplot")
library("viridis")

source("/dfs/Liston_Lab/workspace/caroline/programs/R_plots_themes.R")

## Analysis Cov, Fst, pi for Fvb6-1 only 
mydir = "/dfs/Liston_Lab/scratch/cauretc/2021_sex_chrom_analysis/wgs/pop_stats/13_F_M_no_miss_filt"

#subgenomes
B1 <- c("Fvb1-1", "Fvb2-1", "Fvb3-1", "Fvb4-1", "Fvb5-4", "Fvb6-2", "Fvb7-4")
Bi <- c("Fvb1-2", "Fvb2-4", "Fvb3-2", "Fvb4-4", "Fvb5-3", "Fvb6-3", "Fvb7-3")
B2 <- c("Fvb1-3", "Fvb2-3", "Fvb3-3", "Fvb4-2", "Fvb5-2", "Fvb6-4", "Fvb7-1")
Av <- c("Fvb1-4", "Fvb2-2", "Fvb3-4", "Fvb4-3", "Fvb5-1", "Fvb6-1", "Fvb7-2")

Fvb1 <- c("Fvb1-1", "Fvb1-2", "Fvb1-3", "Fvb1-4")
Fvb2 <- c("Fvb2-1", "Fvb2-2", "Fvb2-3", "Fvb2-4")
Fvb3 <- c("Fvb3-1", "Fvb3-2", "Fvb3-3", "Fvb3-4")
Fvb4 <- c("Fvb4-1", "Fvb4-2", "Fvb4-3", "Fvb4-4")
Fvb5 <- c("Fvb5-1", "Fvb5-2", "Fvb5-3", "Fvb5-4")
Fvb6 <- c("Fvb6-1", "Fvb6-2", "Fvb6-3", "Fvb6-4")
Fvb7 <- c("Fvb7-1", "Fvb7-2", "Fvb7-3", "Fvb7-4")

#Fst
Fstfiles = list.files(path=mydir, pattern="*.windowed.weir.fst", full.names=TRUE)
Fst.10kb = ldply(Fstfiles, read.table, header = TRUE) %>% 
  mutate(CHROM = str_remove_all(CHROM, "_RagTag")) %>% 
  mutate(subgenome = case_when( 
    CHROM %in% unlist(B1) ~ "B1",
    CHROM %in% unlist(Bi) ~ "Bi",
    CHROM %in% unlist(B2) ~ "B2",
    CHROM %in% unlist(Av) ~ "Av",
    TRUE ~ "unknown")) %>%
  mutate(homeologs = case_when( CHROM %in% unlist(Fvb1) ~ "Fchil1",
                                CHROM %in% unlist(Fvb2) ~ "Fchil2",
                                CHROM %in% unlist(Fvb3) ~ "Fchil3",
                                CHROM %in% unlist(Fvb4) ~ "Fchil4",
                                CHROM %in% unlist(Fvb5) ~ "Fchil5",
                                CHROM %in% unlist(Fvb6) ~ "Fchil6",
                                CHROM %in% unlist(Fvb7) ~ "Fchil7",
                                TRUE ~ "unknown"))

Fst_Fc6_1 <- Fst.10kb %>% filter(CHROM=="Fvb6-1") %>% mutate(analysis = "Intersexual Fst") %>% 
  mutate(value = WEIGHTED_FST)%>% select(BIN_START, value, analysis)

#Coverage
all_depth_wind <- fread("/dfs/Liston_Lab/scratch/cauretc/2021_sex_chrom_analysis/wgs/depth_files/no_secondary/ratio_merged_window.txt", sep="\t")
names(all_depth_wind) = gsub(pattern = "#", replacement = "", x = names(all_depth_wind))
cov_sex_chr_wind <- all_depth_wind %>% filter(CHROM == "Fvb6-1_RagTag") %>% mutate(FEMALE_AVG = FEMALE_AVG + 0.001, MALE_AVG = MALE_AVG + 0.001)
cov_Fc6_1 <- cov_sex_chr_wind %>% 
  #mutate(analysis = "Cov F/M") %>%
  mutate(analysis = "Coverage (F:M)") %>%
  mutate(BIN_START = POS_START) %>% 
  #mutate(value = FEMALE_AVG/MALE_AVG)
  mutate(value =log2(FEMALE_AVG)-log2(MALE_AVG))%>% 
  select(BIN_START, value, analysis)

#Merging both analysis
allFc6_1_analysis <- bind_rows(Fst_Fc6_1,cov_Fc6_1)

#Getting some CI (from Fst.R and cov.R)

#Fst - 13FM - no missing filt, weighted
#highCI
#[1] 0.0192476
#> lowCI
#[1] -0.0228739

#Coverage - ratio
#> lowCI 
#[1] 0.8689253
#> highCI 
#[1] 1.155434

lowCI_Fst <- -0.0228739
highCI_Fst <- 0.0192476
lowCI_cov <- -0.2026959
highCI_cov <- 0.1916831

CI_table <- data_frame(low = c(lowCI_cov,lowCI_Fst), high = c(highCI_cov,highCI_Fst), analysis =c("Coverage (F:M)","Intersexual Fst"))

##Add line to highlight the same region plotted in first fig
line_table <- data_frame(start = 33.8, end = 35.3, analysis ="Coverage (F:M)")

#Plotting - each analysis corresponding to a facet
allFc6_1_analysis.plot <- ggplot(allFc6_1_analysis) + 
  scale_x_continuous(name="Fchil6-Av Position (Mbp)") +
  #annotate(data=CI_table, "rect", xmin=-Inf, xmax=Inf, ymin=low, ymax=high, alpha=0.8, fill="grey")+
  geom_rect(data=CI_table, aes(xmin=-Inf, xmax=Inf, ymin=low, max=high), fill="grey", alpha=0.5) +
  geom_rect(data=line_table, aes(xmin=start, xmax=end, ymin=-1.53, max=-1.47), fill="black") +
  #geom_vline(xintercept = 34.2, color="purple", alpha=0.5) +
  #scale_y_continuous(name="Pi.F-Pi.M") +
  geom_line(aes(x=BIN_START/1000000, y= value), colour = "red", size = 0.2)+
  facet_grid(rows=vars(analysis), scales="free") +
  cleanPlot2() +
  theme(axis.title.y=element_blank(),
        strip.text.y = element_text(angle = -90, size = 14))

#ggsave(allFc6_1_analysis.plot,filename = "/dfs/Liston_Lab/workspace/caroline/2021_sex_chrom_analysis/figures/all_10kb_analysis_plot_equalMF.jpeg")

##Autosome 3-1
Fst_Fc3_1 <- Fst.10kb %>% filter(CHROM=="Fvb3-1") %>% mutate(analysis = "Intersexual Fst") %>% 
  mutate(value = WEIGHTED_FST)%>% select(BIN_START, value, analysis)
cov_autosome <- all_depth_wind %>% filter(CHROM == "Fvb3-1_RagTag") %>% 
  mutate(FEMALE_AVG = FEMALE_AVG + 0.001, MALE_AVG = MALE_AVG + 0.001, ratio = FEMALE_AVG/MALE_AVG, log_diff = log2(FEMALE_AVG) - log2(MALE_AVG))
cov_Fc3_1 <- cov_autosome %>% 
  #mutate(analysis = "Cov F/M") %>%
  mutate(analysis = "Coverage (F:M)") %>%
  mutate(BIN_START = POS_START) %>% 
  #mutate(value = FEMALE_AVG/MALE_AVG) %>% 
  mutate(value =log2(FEMALE_AVG)-log2(MALE_AVG))%>% 
  select(BIN_START, value, analysis)
div_Fc3_1 <- divF.10kb %>% filter(CHROM=="Fvb3-1") %>% mutate(analysis = "Pi F") %>% 
  mutate(value = PI) %>% select(BIN_START, value, analysis)
#allFc3_1_analysis <- bind_rows(Fst_Fc3_1,cov_Fc3_1,div_Fc3_1)
allFc3_1_analysis <- bind_rows(Fst_Fc3_1,cov_Fc3_1)

#Merging the analysis from the sex chromosome and autosome
allFc6_1_analysis <- allFc6_1_analysis %>% mutate(chr = "Fchil6-Av")
allFc3_1_analysis <- allFc3_1_analysis %>% mutate(chr = "Fchil3-B1")
all_analysis_SC_autosome <- bind_rows(allFc6_1_analysis,allFc3_1_analysis)

#allFc3_1_analysis.plot <- ggplot(allFc3_1_analysis) + 
#  scale_x_continuous(name="Fchil3-B1 Position (Mbp)") +
#  #geom_rect(data=CI_table, aes(xmin=-Inf, xmax=Inf, ymin=low, max=high), fill="grey", alpha=0.5) +
#  geom_line(aes(x=BIN_START/1000000, y= value), colour = "red", size = 0.2)+
#  facet_grid(rows=vars(analysis), scales="free") +
#  cleanPlot2() +
#  theme(axis.title.y=element_blank(),
#        strip.text.y = element_text(angle = -90, size = 14))

all_analysis_SC_autosome.plot <- ggplot(all_analysis_SC_autosome) + 
  scale_x_continuous(name="Genomic Position (Mbp)") +
  geom_rect(data=CI_table, aes(xmin=-Inf, xmax=Inf, ymin=low, max=high), fill="grey", alpha=0.5) +
  geom_line(aes(x=BIN_START/1000000, y= value, colour = chr), size = 0.2)+
  facet_grid(rows=vars(analysis), cols=vars(chr), scales="free") +
  scale_color_manual(values = c("Fchil6-Av" = "red", "Fchil3-B1" = "#CCBB44"), guide="none")+
  cleanPlot2() +
  theme(axis.title.y=element_blank(),
        strip.text.y = element_text(angle = -90, size = 14))

#ggsave(allFc3_1_analysis.plot,filename = "/dfs/Liston_Lab/workspace/caroline/2021_sex_chrom_analysis/figures/all_10kb_analysis_plot_equalMF_no_CI_only_F_div_autosome.jpeg")
#ggsave(all_analysis_SC_autosome.plot,filename = "/dfs/Liston_Lab/workspace/caroline/2021_sex_chrom_analysis/figures/all_10kb_analysis_plot_equalMF_no_CI_only_F_div_SC_autosome.jpeg")
ggsave(all_analysis_SC_autosome.plot,filename = "/dfs/Liston_Lab/workspace/caroline/2021_sex_chrom_analysis/figures/all_10kb_analysis_plot_equalMF_SC_autosome.jpeg")

#Adding the TEs on the sex chromosome analysis plot
TE_cov_wind <- read.table("/dfs/Liston_Lab/workspace/sebastian/gp33_repeats/20211019_repeat_correction/20211019_rt.scaffs.curated.reorient.wW.fa.out.SexContigs_GP33_10kb_windows.COV.txt", sep ="\t", h = F)
colnames(TE_cov_wind) <- c("chrom", "start", "end", "features_numb", "bp", "wind_size", "proportion")
TE_cov_wind_SC <- TE_cov_wind %>% filter(chrom == "Fvb6-1_RagTag")

TE_grad_plot_SC <- ggplot(TE_cov_wind_SC) +
  geom_tile(aes(y = "SC", x=start/1000000, fill = proportion)) +
  scale_fill_viridis(name = "Repeat cov") +
  #expand_limits(x = c(0,1.4)) +
  #scale_y_discrete(position = "left")+
  scale_y_discrete(position = "right")+
  ylab('Rep') +
  cleanPlot2() +
  scale_x_continuous(position = "top") +
  theme(legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12)) 

legend_TE <- cowplot::get_legend(TE_grad_plot_SC) 

TE_grad_plot_SC <- TE_grad_plot_SC + theme(legend.position = "none", #axis.text.x = element_blank(),  
                                         axis.title.x=element_blank(), axis.ticks.y = element_blank(),
                                         #axis.text = element_text(size = 9),
                                         axis.text.y = element_blank(),
                                         axis.title = element_text(size = 13),
                                         strip.text.y = element_text(angle = -90, size = 14))

first_col = plot_grid(TE_grad_plot_SC, allFc6_1_analysis.plot, rel_heights = c(0.3, 3), ncol = 1, align = "v", axis = "lr")
#ggsave(first_col,filename = "/dfs/Liston_Lab/workspace/caroline/2021_sex_chrom_analysis/figures/all_10kb_analysis_plot_equalMF_SC_TE.jpeg")
ggsave(first_col,filename = "/dfs/Liston_Lab/workspace/caroline/2021_sex_chrom_analysis/figures/all_10kb_analysis_plot_equalMF_SC_TE_line.jpeg")

##After the 1st round of review it was suggested to zoom in at the Fst spike.
#34170001 to 34200001 (30kb)

#TE_cov_wind_SC_subset <- TE_cov_wind_SC %>% filter(start >= 34170000 & start <= 34200000)
#allFc6_1_analysis_subset <- allFc6_1_analysis %>% filter(BIN_START >= 34170001 & BIN_START <= 34200001)
TE_cov_wind_SC_subset <- TE_cov_wind_SC %>% mutate(start = start +1) %>% filter(start >= 34140001 & start <= 34240001)
allFc6_1_analysis_subset <- allFc6_1_analysis %>% filter(BIN_START >= 34140001 & BIN_START <= 34240001)

allFc6_1_analysis_subset.plot <- ggplot(allFc6_1_analysis_subset) + 
  #scale_x_continuous(name="Fchil6-Av Position (Mbp)") +
  #scale_x_continuous(name="Fchil6-Av Position (Mbp)", limits = c(34.13,34.25), breaks = c(34.14, 34.20,34.25)) +
  #scale_x_continuous(name="Fchil6-Av Position (Mbp)", limits = c(34.13,34.25), breaks = c(34.14, 34.19, 34.20, 34.25),
  #                   labels = c("34.14", "SDR", "34.20", "34.25")) +
  scale_x_continuous(name="Fchil6-Av Position (Mbp)", limits = c(34.13,34.25), breaks = c(34.15, 34.19, 34.25),
                     labels = c("34.15", "SDR", "34.25")) +
  #scale_y_continuous(position = "right")+
  #annotate(data=CI_table, "rect", xmin=-Inf, xmax=Inf, ymin=low, ymax=high, alpha=0.8, fill="grey")+
  geom_rect(data=CI_table, aes(xmin=-Inf, xmax=Inf, ymin=low, max=high), fill="grey", alpha=0.5) +
  #geom_rect(data=line_table, aes(xmin=start, xmax=end, ymin=-1.53, max=-1.47), fill="black") +
  #geom_vline(xintercept = 34.2, color="purple", alpha=0.5) +
  #scale_y_continuous(name="Pi.F-Pi.M") +
  geom_line(aes(x=BIN_START/1000000, y= value), colour = "red", size = 0.2)+
  #expand_limits(x = c(34.14,34.25)) +
  #xlim(34.14,34.24) +
  facet_grid(rows=vars(analysis), scales="free") +
  cleanPlot2() +
  theme(axis.title.y=element_blank(),
        strip.text.y = element_text(angle = -90, size = 14),
        #axis.text.x = element_text(color = c("black", "red", "black", "black", "black")),
        #axis.ticks.x = element_line(color = c("black", "red", "black", "black", "black")))
        axis.text.x = element_text(color = c("black", "red", "black")),
        axis.ticks.x = element_line(color = c("black", "red", "black")))

TE_grad_plot_SC_subset <- ggplot(TE_cov_wind_SC_subset) +
  geom_tile(aes(y = "SC", x=start/1000000, fill = proportion)) +
  scale_fill_viridis(name = "Repeat cov") +
  #expand_limits(x = c(34.14,34.24)) +
  #xlim(34.14,34.24) +
  #scale_y_discrete(position = "left")+
  scale_y_discrete(position = "right")+
  ylab('Rep') +
  cleanPlot2() +
  #scale_x_continuous(position = "top") +
  #scale_x_continuous(position = "top", limits = c(34.13,34.25), breaks = c(34.14, 34.20,34.25)) +
  scale_x_continuous(position = "top", limits = c(34.13,34.25), breaks = c(34.15, 34.19, 34.25),
                     labels = c("34.15", "SDR", "34.25")) +
  #scale_x_continuous(position = "top", limits = c(34.13,34.25), breaks = c(34.14, 34.19, 34.20, 34.25),
                     #labels = c("34.14", "SDR", "34.20", "34.25")) +
  theme(legend.key.size = unit(0.5, 'cm'),
        #axis.text.x = element_text(color = c("black", "red", "black", "black", "black")),
        #axis.ticks.x = element_line(color = c("black", "red", "black", "black", "black")),
        axis.text.x = element_text(color = c("black", "red", "black", "black")),
        axis.ticks.x = element_line(color = c("black", "red", "black", "black")),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12)) 

#legend_TE <- cowplot::get_legend(TE_grad_plot_SC) 

TE_grad_plot_SC_subset <- TE_grad_plot_SC_subset + theme(legend.position = "none", #axis.text.x = element_blank(),  
                                           axis.title.x=element_blank(), axis.ticks.y = element_blank(),
                                           #axis.text = element_text(size = 9),
                                           axis.text.y = element_blank(),
                                           axis.title = element_text(size = 13),
                                           strip.text.y = element_text(angle = -90, size = 14))

#first_col = plot_grid(TE_grad_plot_SC_subset, allFc6_1_analysis_subset.plot, rel_heights = c(0.3, 3), ncol = 1, align = "v", axis = "lr")
#ggsave(first_col,filename = "/dfs/Liston_Lab/workspace/caroline/2021_sex_chrom_analysis/figures/all_10kb_analysis_plot_equalMF_SC_TE_zoom_Fstspike.jpeg")
#ggsave(first_col,filename = "/dfs/Liston_Lab/workspace/caroline/2021_sex_chrom_analysis/figures/all_10kb_analysis_plot_equalMF_SC_TE_zoom_Fstspike_test.jpeg")

#combine both: all sex chromosome and zoomed as separate columns
allFc6_1_analysis.plot <- ggplot(allFc6_1_analysis) + 
  scale_x_continuous(name="Fchil6-Av Position (Mbp)") +
  geom_rect(data=CI_table, aes(xmin=-Inf, xmax=Inf, ymin=low, max=high), fill="grey", alpha=0.5) +
  geom_rect(data=line_table, aes(xmin=start, xmax=end, ymin=-1.53, max=-1.47), fill="black") +
  geom_line(aes(x=BIN_START/1000000, y= value), colour = "red", size = 0.2)+
  facet_grid(rows=vars(analysis), scales="free") +
  cleanPlot2() +
  theme(axis.title.y=element_blank(),
        #strip.text.y = element_text(angle = -90, size = 14)
        strip.text.y = element_blank())

TE_grad_plot_SC <- TE_grad_plot_SC + theme(legend.position = "none", #axis.text.x = element_blank(),  
                                           axis.title.x=element_blank(), axis.ticks.y = element_blank(),
                                           #axis.text = element_text(size = 9),
                                           axis.text.y = element_blank(),
                                           #axis.title = element_text(size = 13),
                                           axis.title = element_blank(),
                                           #strip.text.y = element_text(angle = -90, size = 14)
                                           strip.text.y = element_blank())
first_col = plot_grid(TE_grad_plot_SC, allFc6_1_analysis.plot, rel_heights = c(0.3, 3), ncol = 1, align = "v", axis = "lr")
second_col = plot_grid(TE_grad_plot_SC_subset, allFc6_1_analysis_subset.plot, rel_heights = c(0.3, 3), ncol = 1, align = "v", axis = "lr")
perfect = plot_grid(first_col, second_col, ncol = 2, rel_widths = c(1, 1)) + theme(plot.background = element_rect(fill = "white"))
ggsave(perfect,filename = "/dfs/Liston_Lab/workspace/caroline/2021_sex_chrom_analysis/figures/all_10kb_analysis_plot_equalMF_SC_TE_zoom_Fstspike_test_both.jpeg")
