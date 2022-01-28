# Creating Figure Fst figures: Genome wide with all populations combined and subset (sex chromosome) with population separated
# This was used to produce the Fig S8, S9, S10 of Cauret et al. 2022
# Chromosome-scale assembly with a phased sex-determining region resolves features of early Z and W chromosome differentiation in a wild octoploid strawberry

# Loading libraries
library("dplyr")
library("plyr")
library("ggplot2")
library("stringr")
library("grid")

source("/dfs/Liston_Lab/workspace/caroline/programs/R_plots_themes.R")

#Genome wide Fst
#mydir = "/dfs/Liston_Lab/scratch/cauretc/2021_sex_chrom_analysis/wgs/pop_stats" #when all individuals were used
mydir = "/dfs/Liston_Lab/scratch/cauretc/2021_sex_chrom_analysis/wgs/pop_stats/13_F_M_no_miss_filt"
Fstfiles = list.files(path=mydir, pattern="*.windowed.weir.fst", full.names=TRUE)

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


Fst.10kb = ldply(Fstfiles, read.table, header = TRUE) %>% 
  mutate(CHROM = str_remove_all(CHROM, "_RagTag")) %>% 
  mutate(subgenome = case_when( #CHROM %in% unlist(vir) ~ "viridis",
    #CHROM %in% unlist(iin) ~ "iinumae",
    #CHROM %in% unlist(nip) ~ "nipponica",
    #CHROM %in% unlist(ves) ~ "vesca",
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
Fst_Fc3_1 <- Fst.10kb%>% filter(CHROM=="Fvb3-1")
myfunction <- function(i){
  Info <- sample(i,1,replace=FALSE)
  return(Info)
}

sampling_number <-1000
my.perm <- c()

for(i in 1:sampling_number){ my.perm[i] <- myfunction(Fst_Fc3_1$WEIGHTED_FST) }
#for(i in 1:sampling_number){ my.perm[i] <- myfunction(Fst_Fc3_1$MEAN_FST) }
sorted.perm <- sort(my.perm)
lowCI <- sorted.perm[25]
highCI <- sorted.perm[975]

#13FM - no missing filt, weighted
#highCI
#[1] 0.0192476
#> lowCI
#[1] -0.0228739

Fst.10kb.plot_CI <- ggplot(Fst.10kb) + 
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=lowCI, ymax=highCI, alpha=0.8, fill="grey") +
  geom_line(aes(x=BIN_START/1000000, y=WEIGHTED_FST, color = subgenome))+
  #geom_line(aes(x=BIN_START/1000000, y=MEAN_FST, color = subgenome))+
  scale_x_continuous(name="Genomic Position (Mbp)") +
  scale_color_manual(values = c("Av" = "red", "B1" = "#CCBB44", 
                                "B2" = "#228833", "Bi" = "#4477AA"), guide="none")+
  facet_grid(rows=vars(homeologs), cols = vars(subgenome)) +
  coord_cartesian(ylim = c(-0.06, 0.22))+
  cleanPlot2() 

#ggsave(Fst.10kb.plot_CI,filename = "/dfs/Liston_Lab/scratch/cauretc/2021_sex_chrom_analysis/wgs/pop_stats/Fst_10kb_plot_subg_CI.jpeg")
ggsave(Fst.10kb.plot_CI,filename = "/dfs/Liston_Lab/scratch/cauretc/2021_sex_chrom_analysis/wgs/pop_stats/13_F_M_no_miss_filt/Fst_mean_13_F_M_no_miss_filt_10kb_plot_subg_CI.jpeg")

#Fst per population
Fst_all <- read.table("/dfs/Liston_Lab/scratch/cauretc/2021_sex_chrom_analysis/wgs/pop_stats/13_F_M_no_miss_filt/chil_wgs_Liston_UCD_all_ind_remapped_Fvb6-1_RagTag_filtered_ind_filtered_qual_13_F_VS_13_M.windowed.weir.fst",
                     sep ="\t", h = T) %>% mutate(Population = "All")

Fst_OR <- read.table("/dfs/Liston_Lab/scratch/cauretc/2021_sex_chrom_analysis/wgs/pop_stats/by_pop_OR_CA/Sex_chr_OR_F_VS_M.windowed.weir.fst", 
                     sep ="\t", h = T) %>% mutate(Population = "Oregon")

Fst_CA <- read.table("/dfs/Liston_Lab/scratch/cauretc/2021_sex_chrom_analysis/wgs/pop_stats/by_pop_OR_CA/Sex_chr_CA_F_VS_M.windowed.weir.fst", 
                     sep ="\t", h = T) %>% mutate(Population = "California")

Fst_all_pop <- bind_rows(Fst_all,Fst_OR,Fst_CA)


Fst_all_pop.plot <- ggplot(Fst_all_pop) + 
  scale_x_continuous(name="Fchil6-Av Position (Mbp)") +
  #annotate(data=CI_table, "rect", xmin=-Inf, xmax=Inf, ymin=low, ymax=high, alpha=0.8, fill="grey")+
  #geom_rect(data=CI_table, aes(xmin=-Inf, xmax=Inf, ymin=low, max=high), fill="grey", alpha=0.5) +
  #geom_vline(xintercept = 34.2, color="purple", alpha=0.5) +
  #scale_y_continuous(name="Pi.F-Pi.M") +
  geom_line(aes(x=BIN_START/1000000, y=WEIGHTED_FST), colour = "red", size = 0.2)+
  facet_grid(rows=vars(Population), scales="free") +
  cleanPlot2() 

ggsave(Fst_all_pop.plot,filename = "/dfs/Liston_Lab/workspace/caroline/2021_sex_chrom_analysis/figures/Fst_sex_chr_all_OR_CA.jpeg")



