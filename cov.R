# Creating Figure combining read coverage, nucmer alignment of ZW sequence and TE distribution
# This was used to produce the main figure 3 of Cauret et al. 2022
# Chromosome-scale assembly with a phased sex-determining region resolves features of early Z and W chromosome differentiation in a wild octoploid strawberry

# Loading libraries
library("dplyr")
require("data.table")
library("tidyr")
library("zoo")
library("ggplot2")
library("grid")
library("cowplot")
library("viridis")

source("/dfs/Liston_Lab/workspace/caroline/programs/R_plots_themes.R")
all_depth_wind <- fread("/dfs/Liston_Lab/scratch/cauretc/2021_sex_chrom_analysis/wgs/depth_files/no_secondary/ratio_merged_window.txt", sep="\t")
names(all_depth_wind) = gsub(pattern = "#", replacement = "", x = names(all_depth_wind))

cov_W_wind <- all_depth_wind %>% filter(CHROM == "atg000165l") %>% mutate(FEMALE_AVG = FEMALE_AVG + 0.001, MALE_AVG = MALE_AVG + 0.001)
cov_sex_chr_wind <- all_depth_wind %>% filter(CHROM == "Fvb6-1_RagTag") %>% mutate(FEMALE_AVG = FEMALE_AVG + 0.001, MALE_AVG = MALE_AVG + 0.001)

##Confidence intervals
cov_autosome <- all_depth_wind %>% filter(CHROM == "Fvb3-1_RagTag") %>% 
  mutate(FEMALE_AVG = FEMALE_AVG + 0.001, MALE_AVG = MALE_AVG + 0.001, ratio = FEMALE_AVG/MALE_AVG, log_diff = log2(FEMALE_AVG) - log2(MALE_AVG))
  

myfunction <- function(i){
  Info <- sample(i,1,replace=FALSE)
  return(Info)
}

sampling_number <-1000
my.perm <- c()
#for(i in 1:sampling_number){ my.perm[i] <- myfunction(cov_autosome$ratio) }
for(i in 1:sampling_number){ my.perm[i] <- myfunction(cov_autosome$log_diff) }
sorted.perm <- sort(my.perm)
lowCI <- sorted.perm[25]
highCI <- sorted.perm[975]
#ratio
#> lowCI 
#[1] 0.8689253
#> highCI 
#[1] 1.155434
#log2
#> lowCI 
#[1] -0.2026959
#> highCI 
#[1] 0.1916831

lowCI <- -0.2026959
highCI <- 0.1916831

##nucmer ZW
readDelta <- function(deltafile){
  lines = scan(deltafile, 'a', sep='\n', quiet=TRUE)
  lines = lines[-1]
  lines.l = strsplit(lines, ' ')
  lines.len = lapply(lines.l, length) %>% as.numeric
  lines.l = lines.l[lines.len != 1]
  lines.len = lines.len[lines.len != 1]
  head.pos = which(lines.len == 4)
  head.id = rep(head.pos, c(head.pos[-1], length(lines.l)+1)-head.pos)
  mat = matrix(as.numeric(unlist(lines.l[lines.len==7])), 7)
  res = as.data.frame(t(mat[1:5,]))
  colnames(res) = c('rs','re','qs','qe','error')
  res$qid = unlist(lapply(lines.l[head.id[lines.len==7]], '[', 2))
  res$rid = unlist(lapply(lines.l[head.id[lines.len==7]], '[', 1)) %>% gsub('^>', '', .)
  res$strand = ifelse(res$qe-res$qs > 0, '+', '-')
  res
}

interval <- data.frame(chrom="atg000165l", start=298666, end=334844)
interval_Z <- data.frame(chrom="Fvb6-1", start=34170001, end=34190001)

mumgp = readDelta("/dfs/Liston_Lab/workspace/caroline/2021_sex_chrom_analysis/nucmer_W_Z/Fvb61_atgW_filter.delta")

mumplot <- ggplot(mumgp, aes(x=qs/1000000, xend=qe/1000000, y=rs/1000000, yend=re/1000000, colour=strand)) + geom_segment(size=1.5) +
  theme_bw() + cleanPlot2() + ylim(33.8,35.3) + theme(legend.position = "none",
                                                      axis.text.y = element_text(),
                                                      axis.text = element_text(size = 9),
                                                      axis.title = element_text(size = 13)
  ) +
  xlab('Position on W haplotig (Mb)') + ylab('Position on Fchil6-Av (Mb)') + 
  scale_x_continuous(position = "top") +
  expand_limits(x = c(0,1.4)) + scale_color_manual(values = c("blue", "red")) +
  geom_rect(data=interval, inherit.aes=FALSE, aes(xmin=start/1000000, xmax=end/1000000, ymin=-Inf,
                                                  ymax=Inf), color="transparent", fill="purple", alpha=0.5) +
  geom_rect(data=interval_Z, inherit.aes=FALSE, aes(ymin=start/1000000, ymax=end/1000000, xmin=-Inf,
                                                  xmax=Inf), color="transparent", fill="violet", alpha=0.5)

depth.average.plot.W <- ggplot(cov_W_wind) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=lowCI, ymax=highCI, alpha=0.8, fill="grey") +
  geom_rect(data=interval, inherit.aes=FALSE, aes(xmin=start/1000000, xmax=end/1000000, ymin=-Inf,
                                                  ymax=Inf), color="transparent", fill="purple", alpha=0.5) +
  #geom_line(aes(x=POS_START/1000000, y=FEMALE_AVG/MALE_AVG)) +
  geom_line(aes(x=POS_START/1000000, y=log2(FEMALE_AVG)-log2(MALE_AVG))) +
  #theme_bw() + 
  cleanPlot2()+ theme(legend.position = "none",axis.text.y = element_text(),  axis.title.x=element_blank(),
                      axis.text = element_text(size = 9),
                      axis.title = element_text(size = 13)) +
  xlab('W contig') + 
  #ylab('Cov F/M') + 
  #ylab('Log2(F:M)') +
  ylab('Cov(F:M)') +
  expand_limits(x = c(0,1.4)) 

depth.average.plot.SC <- ggplot(cov_sex_chr_wind)+
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=lowCI, ymax=highCI, alpha=0.8, fill="grey") +
  geom_rect(data=interval_Z, inherit.aes=FALSE, aes(xmin=start/1000000, xmax=end/1000000, ymin=-Inf,
                                                    ymax=Inf), color="transparent", fill="violet", alpha=0.5) +
  #geom_line(aes(x=POS_START/1000000, y=FEMALE_AVG/MALE_AVG)) + 
  geom_line(aes(x=POS_START/1000000, y=log2(FEMALE_AVG)-log2(MALE_AVG))) +
  cleanPlot2() +
  theme(#legend.position = "none",
    axis.title.y=element_blank(),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 13)
  )  + xlab('Fchil6-Av') + scale_y_continuous(position = "right") +
  #ylab('Cov F/M')+
  #ylab('Log2(F:M)') +
  ylab('Cov(F:M)') +
  scale_x_continuous(position = "top", limits=c(33.8,35.3)) 

depth.average.plot.SC <- depth.average.plot.SC + coord_flip() 
#changing margin
#mumplot = mumplot + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
#depth.average.plot.W = depth.average.plot.W + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
#depth.average.plot.SC = depth.average.plot.SC + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

#first_col = plot_grid(mumplot, depth.average.plot.W, ncol = 1, rel_heights = c(3, 1), align = "hv")
#second_col = plot_grid(depth.average.plot.SC,NULL, ncol = 1, rel_heights = c(3, 1), align = "hv")
#perfect = plot_grid(first_col, second_col, ncol = 2, rel_widths = c(3, 1), nrow = 1, align = "hv") + theme(plot.background = element_rect(fill = "white"))
first_col = plot_grid(mumplot, depth.average.plot.W, rel_heights = c(3, 1), ncol = 1, align = "v")
second_col = plot_grid(depth.average.plot.SC,NULL, rel_heights = c(3, 1), ncol = 1)
perfect = plot_grid(first_col, second_col, ncol = 2, rel_widths = c(3, 1)) + theme(plot.background = element_rect(fill = "white"))


#perfect = plot_grid(mumplot, depth.average.plot.SC, depth.average.plot.W, ncol = 2, rel_widths = c(3, 1), nrow = 2)
#perfect = plot_grid(mumplot, depth.average.plot.SC, depth.average.plot.W, ncol = 2, rel_heights = c(3, 1), rel_widths = c(3, 1), nrow = 2, align = "hv") + theme(plot.background = element_rect(fill = "white"))

#ggsave(perfect,filename = "/dfs/Liston_Lab/workspace/caroline/2021_sex_chrom_analysis/figures/mumplot_ZW_cov_combined_no_sec.jpeg")
ggsave(perfect,filename = "/dfs/Liston_Lab/workspace/caroline/2021_sex_chrom_analysis/figures/mumplot_ZW_cov_combined_no_sec_log.jpeg")

#Add TE gradient
TE_cov_wind <- read.table("/dfs/Liston_Lab/workspace/sebastian/gp33_repeats/20211019_repeat_correction/20211019_rt.scaffs.curated.reorient.wW.fa.out.SexContigs_GP33_10kb_windows.COV.txt", sep ="\t", h = F)
colnames(TE_cov_wind) <- c("chrom", "start", "end", "features_numb", "bp", "wind_size", "proportion")

TE_cov_wind_Z <- TE_cov_wind %>% filter(chrom == "Fvb6-1_RagTag")
TE_cov_wind_W <- TE_cov_wind %>% filter(chrom == "atg000165l")

TE_grad_plot_Z <- ggplot(TE_cov_wind_Z) +
  #geom_rect(aes(ymin = 1, ymax = 1, xmin=min(start/1000000), xmax= max(start/1000000), fill = proportion)) +
  #geom_tile(aes(y = 1, x=start/1000000, fill = proportion)) +
  #geom_tile(aes(y = start/1000000, x=1, fill = proportion)) +
  geom_tile(aes(y = start/1000000, x="Z", fill = proportion)) +
  #scale_y_continuous(limits=c(0.5,1.5),breaks=1)+
  #scale_x_continuous(limits=c(0.5,1.5),breaks=1, position = "top")+
  #scale_fill_gradient(low = 'blue', high = 'red') +
  scale_fill_viridis() +
  #scale_fill_gradient2(low = 'blue', high = 'red', midpoint = 0.5) +
  #theme_bw()
  #ylab('TE') +
  #expand_limits(x=c(33.8,35.3))+
  #xlab('TE') +
  xlab('Rep') +
  ylim(y=c(33.8,35.3))+
  scale_x_discrete(position = "top")+
  #expand_limits(y=c(33.8,35.3))+
  cleanPlot2() +
  #theme(legend.position = "none")
  theme(legend.position = "none", axis.text.y = element_blank(), 
        axis.title.y=element_blank(), axis.ticks.y = element_blank(),
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 13))

#legend_TE_1 <- cowplot::get_legend(TE_grad_plot_Z) 
#ggsave(legend_TE_1,filename = "/dfs/Liston_Lab/workspace/caroline/2021_sex_chrom_analysis/figures/TE_grad__Z_legend_check.jpeg")

#TE_grad_plot_Z <- TE_grad_plot_Z + coord_flip() 

TE_grad_plot_W <- ggplot(TE_cov_wind_W) +
  #geom_tile(aes(y = 1, x=start/1000000, fill = proportion)) +
  geom_tile(aes(y = "W", x=start/1000000, fill = proportion)) +
  #geom_rect(aes(ymin = 1, ymax = 1, xmin=min(start/1000000), xmax= max(start/1000000), fill = proportion)) +
  #scale_y_continuous(limits=c(0.5,1.5),breaks=1)+
  #scale_fill_gradient(low = 'blue', high = 'red') +
  scale_fill_viridis(name = "Repeat cov") +
  #scale_fill_gradient2(low = 'blue', high = 'red', midpoint = 0.5) +
  expand_limits(x = c(0,1.4)) +
  scale_y_discrete(position = "left")+
  #ylab('TE') +
  ylab('Rep') +
  cleanPlot2() +
  theme(legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12))
  #theme(legend.position = "none", axis.text.x = element_blank(),  
  #      axis.title.x=element_blank(), axis.ticks.x = element_blank(),
  #      axis.text = element_text(size = 9),
  #      axis.title = element_text(size = 13))

legend_TE <- cowplot::get_legend(TE_grad_plot_W) 

TE_grad_plot_W <- TE_grad_plot_W + theme(legend.position = "none", axis.text.x = element_blank(),  
                                         axis.title.x=element_blank(), axis.ticks.x = element_blank(),
                                         axis.text = element_text(size = 9),
                                         axis.title = element_text(size = 13))

#ggsave(TE_grad_plot_W,filename = "/dfs/Liston_Lab/workspace/caroline/2021_sex_chrom_analysis/figures/TE_grad_try.jpeg")

first_col = plot_grid(mumplot, TE_grad_plot_W, depth.average.plot.W, rel_heights = c(3, 0.2, 1), ncol = 1, align = "v")
second_col = plot_grid(TE_grad_plot_Z,NULL,NULL, rel_heights = c(3, 0.2, 1), ncol = 1)
third_col = plot_grid(depth.average.plot.SC,NULL,legend_TE, rel_heights = c(3, 0.2, 1), ncol = 1)
perfect = plot_grid(first_col, second_col, third_col, ncol = 3, rel_widths = c(3,0.21, 1)) + theme(plot.background = element_rect(fill = "white"))
ggsave(perfect,filename = "/dfs/Liston_Lab/workspace/caroline/2021_sex_chrom_analysis/figures/TE_grad_try.jpeg")

#first_col = plot_grid(mumplot, TE_grad_plot_W, depth.average.plot.W, rel_heights = c(3, 0.2, 1), ncol = 1, align = "v")
#second_col = plot_grid(depth.average.plot.SC,NULL,NULL, rel_heights = c(3, 0.2, 1), ncol = 1)
#perfect = plot_grid(first_col, second_col, ncol = 2, rel_widths = c(3, 1)) + theme(plot.background = element_rect(fill = "white"))
#ggsave(perfect,filename = "/dfs/Liston_Lab/workspace/caroline/2021_sex_chrom_analysis/figures/TE_grad_try.jpeg")

#Z and W spikes
cov_sex_chr_wind_subset <- cov_sex_chr_wind %>% filter(POS_START/1000000 >33.8) %>% mutate(logdiff = log2(FEMALE_AVG)-log2(MALE_AVG))
cov_sex_chr_wind_subset%>%filter(logdiff == min(cov_sex_chr_wind_subset$logdiff))
#Fvb6-1_RagTag  34180001 34190000  0.3035271 0.6447965 -1.087018
cov_sex_chr_wind_subset%>%filter(POS_START >= 34160001 & POS_START <= 34200001)
#CHROM POS_START  POS_END FEMALE_AVG  MALE_AVG     logdiff
#1: Fvb6-1_RagTag  34160001 34170000  0.5432443 0.6385070 -0.23310128
#2: Fvb6-1_RagTag  34170001 34180000  0.4502502 0.8794204 -0.96582616
#3: Fvb6-1_RagTag  34180001 34190000  0.3035271 0.6447965 -1.08701831
#4: Fvb6-1_RagTag  34190001 34200000  0.3668872 0.6762393 -0.88219719
#5: Fvb6-1_RagTag  34200001 34210000  1.5254700 1.4673810  0.05601031

cov_sex_chr_wind_subset <- cov_sex_chr_wind %>% filter(POS_START/1000000 <10) %>% mutate(logdiff = log2(FEMALE_AVG)-log2(MALE_AVG))
cov_sex_chr_wind_subset%>%filter(logdiff == min(cov_sex_chr_wind_subset$logdiff))
#1: Fvb6-1_RagTag    400001  410000  0.1259919 0.3437212 -1.447907
cov_sex_chr_wind_subset%>%filter(POS_START >= 200001 & POS_START <= 700000)

cov_sex_chr_wind <- cov_sex_chr_wind %>% mutate(logdiff = log2(FEMALE_AVG)-log2(MALE_AVG))
cov_sex_chr_wind%>%filter(logdiff == max(cov_sex_chr_wind$logdiff))
#1: Fvb6-1_RagTag  33250001 33260000  0.8617342 0.4373784 0.978361
cov_sex_chr_wind%>%filter(POS_START >= 33230001 & POS_START <= 33270001)
#CHROM POS_START  POS_END FEMALE_AVG  MALE_AVG    logdiff
#1: Fvb6-1_RagTag  33230001 33240000  1.4019561 1.2036429  0.2200337
#2: Fvb6-1_RagTag  33240001 33250000  1.7345344 1.4060333  0.3029177
#3: Fvb6-1_RagTag  33250001 33260000  0.8617342 0.4373784  0.9783610
#4: Fvb6-1_RagTag  33260001 33270000  0.7468994 0.6176946  0.2740202
#5: Fvb6-1_RagTag  33270001 33280000  0.2551614 0.3066997 -0.2654168

cov_W_wind <- cov_W_wind %>% mutate(logdiff = log2(FEMALE_AVG)-log2(MALE_AVG))
cov_W_wind%>%filter(logdiff == max(cov_W_wind$logdiff))
#        CHROM POS_START POS_END FEMALE_AVG  MALE_AVG  logdiff
#1: atg000165l    320001  330000  0.8146771 0.1982684 2.038773
cov_W_wind%>%filter(POS_START >= 270001 & POS_START <= 340001)
#CHROM POS_START POS_END FEMALE_AVG  MALE_AVG    logdiff
#1: atg000165l    270001  280000  0.5951924 0.5753647 0.04887935
#2: atg000165l    280001  290000  0.7318677 0.5794336 0.33693944
#3: atg000165l    290001  300000  0.6612606 0.3286967 1.00846196
#4: atg000165l    300001  310000  1.1219109 0.6678364 0.74839149
#5: atg000165l    310001  320000  1.1450028 0.7251480 0.65900365
#6: atg000165l    320001  330000  0.8146771 0.1982684 2.03877332
#7: atg000165l    330001  340000  0.8206631 0.2067223 1.98909620
#8: atg000165l    340001  350000  0.3901607 0.3458712 0.17383349