# Creating dotplots to compare camarosa and GP33 assembly
# This was used to produce the supplementary figs S1, S2, S3 and main figure 2 of Cauret et al. 2022
# Chromosome-scale assembly with a phased sex-determining region resolves features of early Z and W chromosome differentiation in a wild octoploid strawberry

# Loading libraries
library(pafr)
library(ggplot2)
library(grid)
library(dplyr)
library(gridExtra)
library(stringr)
library(cowplot)
library(patchwork)
library(ggpubr)

# Note: I had to change the pafr program because od the font size which was causing issues
source("/home/caroline/Documents/Caroline_postdoc_OSU/scripts/2021_phil/dotplot_modif_pafr.R")

# Working directory
setwd("/home/caroline/Documents/Caroline_postdoc_OSU/pacbio_chiloensis/minimap2_all/ragtag_round2/manual_curation")

# Alignment between final GP33 assembly and diploid vesca genome
# Renaming the chromosome of chiloensis that was based on camarosa in a more meaningful way, i.e. based on subgenomes
ali_ragtag_curated_scaff_ves <- read_paf("GP33_ragtag_manual_curation_against_F_vesca.paf")%>% 
  mutate(
    qname = str_replace_all( qname, pattern = c("Fvb1-1_RagTag" = "Fchil1-B1", "Fvb1-2_RagTag" = "Fchil1-Bi", "Fvb1-3_RagTag" = "Fchil1-B2", "Fvb1-4_RagTag" ="Fchil1-Av",
                    "Fvb2-1_RagTag" = "Fchil2-B1", "Fvb2-2_RagTag" = "Fchil2-Av", "Fvb2-3_RagTag" = "Fchil2-B2", "Fvb2-4_RagTag" ="Fchil2-Bi",
                    "Fvb3-1_RagTag" = "Fchil3-B1", "Fvb3-2_RagTag" = "Fchil3-Bi", "Fvb3-3_RagTag" = "Fchil3-B2", "Fvb3-4_RagTag" ="Fchil3-Av",
                    "Fvb4-1_RagTag" = "Fchil4-B1", "Fvb4-2_RagTag" = "Fchil4-B2", "Fvb4-3_RagTag" = "Fchil4-Av", "Fvb4-4_RagTag" ="Fchil4-Bi",
                    "Fvb5-1_RagTag" = "Fchil5-Av", "Fvb5-2_RagTag" = "Fchil5-B2", "Fvb5-3_RagTag" = "Fchil5-Bi", "Fvb5-4_RagTag" ="Fchil5-B1",
                    "Fvb6-1_RagTag" = "Fchil6-Av", "Fvb6-2_RagTag" = "Fchil6-B1", "Fvb6-3_RagTag" = "Fchil6-Bi", "Fvb6-4_RagTag" ="Fchil6-B2",
                    "Fvb7-1_RagTag" = "Fchil7-B2", "Fvb7-2_RagTag" = "Fchil7-Av", "Fvb7-3_RagTag" = "Fchil7-Bi", "Fvb7-4_RagTag" ="Fchil7-B1")))


# Ordering the chromosomes
to_keep <- list(
  c("Fchil1-Av", "Fchil1-B1", "Fchil1-Bi", "Fchil1-B2",
    "Fchil2-Av", "Fchil2-B1", "Fchil2-Bi", "Fchil2-B2",
    "Fchil3-Av", "Fchil3-B1", "Fchil3-Bi", "Fchil3-B2",
    "Fchil4-Av", "Fchil4-B1", "Fchil4-Bi", "Fchil4-B2",
    "Fchil5-Av", "Fchil5-B1", "Fchil5-Bi", "Fchil5-B2",
    "Fchil6-Av", "Fchil6-B1", "Fchil6-Bi", "Fchil6-B2",
    "Fchil7-Av", "Fchil7-B1", "Fchil7-Bi", "Fchil7-B2"),
  c("Fvb1", "Fvb2", "Fvb3", "Fvb4", "Fvb5", "Fvb6", "Fvb7")
)

# Alignment between final GP33 assembly and diploid vesca
dot_GP33_vesca <- dotplot1(ali_ragtag_curated_scaff_ves, label_seqs=TRUE, order_by="provided", ordering=to_keep,
                          xlab = "F.chiloensis", ylab = "F.vesca")  +
    theme_bw()

dot_GP33_vesca

ggsave(dot_GP33_vesca,filename = "/home/caroline/Documents/Caroline_postdoc_OSU/pacbio_chiloensis/2021_Phil/FigS1_dotplot_GP33_vesca.jpeg",
       width=20)

# Alignment between final GP33 assembly and octoploid ananassa camarosa genome
## Renaming chromosomes in a more meaningful way, i.e. based on subgenomes
ali_ragtag2round_scaff_cam <- read_paf("GP33_ragtag_manual_curation_against_F_ana_Camarosa.paf")%>% 
  mutate(
    qname = str_replace_all( qname, pattern = c("Fvb1-1_RagTag" = "Fchil1-B1", "Fvb1-2_RagTag" = "Fchil1-Bi", "Fvb1-3_RagTag" = "Fchil1-B2", "Fvb1-4_RagTag" ="Fchil1-Av",
                                                "Fvb2-1_RagTag" = "Fchil2-B1", "Fvb2-2_RagTag" = "Fchil2-Av", "Fvb2-3_RagTag" = "Fchil2-B2", "Fvb2-4_RagTag" ="Fchil2-Bi",
                                                "Fvb3-1_RagTag" = "Fchil3-B1", "Fvb3-2_RagTag" = "Fchil3-Bi", "Fvb3-3_RagTag" = "Fchil3-B2", "Fvb3-4_RagTag" ="Fchil3-Av",
                                                "Fvb4-1_RagTag" = "Fchil4-B1", "Fvb4-2_RagTag" = "Fchil4-B2", "Fvb4-3_RagTag" = "Fchil4-Av", "Fvb4-4_RagTag" ="Fchil4-Bi",
                                                "Fvb5-1_RagTag" = "Fchil5-Av", "Fvb5-2_RagTag" = "Fchil5-B2", "Fvb5-3_RagTag" = "Fchil5-Bi", "Fvb5-4_RagTag" ="Fchil5-B1",
                                                "Fvb6-1_RagTag" = "Fchil6-Av", "Fvb6-2_RagTag" = "Fchil6-B1", "Fvb6-3_RagTag" = "Fchil6-Bi", "Fvb6-4_RagTag" ="Fchil6-B2",
                                                "Fvb7-1_RagTag" = "Fchil7-B2", "Fvb7-2_RagTag" = "Fchil7-Av", "Fvb7-3_RagTag" = "Fchil7-Bi", "Fvb7-4_RagTag" ="Fchil7-B1")))

# Classify by subgenomes
B1 <- c("Fvb1-1", "Fvb2-1", "Fvb3-1", "Fvb4-1", "Fvb5-4", "Fvb6-2", "Fvb7-4")
Bi <- c("Fvb1-2", "Fvb2-4", "Fvb3-2", "Fvb4-4", "Fvb5-3", "Fvb6-3", "Fvb7-3")
B2 <- c("Fvb1-3", "Fvb2-3", "Fvb3-3", "Fvb4-2", "Fvb5-2", "Fvb6-4", "Fvb7-1")
Av <- c("Fvb1-4", "Fvb2-2", "Fvb3-4", "Fvb4-3", "Fvb5-1", "Fvb6-1", "Fvb7-2")


to_keep <- list(
  c("Fchil1-Av", "Fchil1-B1", "Fchil1-Bi", "Fchil1-B2",
    "Fchil2-Av", "Fchil2-B1", "Fchil2-Bi", "Fchil2-B2",
    "Fchil3-Av", "Fchil3-B1", "Fchil3-Bi", "Fchil3-B2",
    "Fchil4-Av", "Fchil4-B1", "Fchil4-Bi", "Fchil4-B2",
    "Fchil5-Av", "Fchil5-B1", "Fchil5-Bi", "Fchil5-B2",
    "Fchil6-Av", "Fchil6-B1", "Fchil6-Bi", "Fchil6-B2",
    "Fchil7-Av", "Fchil7-B1", "Fchil7-Bi", "Fchil7-B2"),
  c("Fvb1-4", "Fvb1-1", "Fvb1-2", "Fvb1-3",
    "Fvb2-2", "Fvb2-1", "Fvb2-4", "Fvb2-3",
    "Fvb3-4", "Fvb3-1", "Fvb3-2", "Fvb3-3",
    "Fvb4-3", "Fvb4-1", "Fvb4-4", "Fvb4-2",
    "Fvb5-1", "Fvb5-4", "Fvb5-3", "Fvb5-2",
    "Fvb6-1", "Fvb6-2", "Fvb6-3", "Fvb6-4",
    "Fvb7-2", "Fvb7-4", "Fvb7-3", "Fvb7-1")
)
dot_GP33_camarosa <- dotplot1(ali_ragtag2round_scaff_cam, label_seqs=TRUE, order_by="provided", ordering=to_keep,
                              xlab = "F.chiloensis", ylab = "F.x camarosa")  +
  theme_bw()
dot_GP33_camarosa
ggsave(dot_GP33_camarosa,filename = "/home/caroline/Documents/Caroline_postdoc_OSU/pacbio_chiloensis/2021_Phil/FigS2_dotplot_GP33_camarosa.jpeg",
       width=20)


# Alignment between final GP33 assembly and octoploid ananassa Royal Royce genome
##RR path: /dfs/Liston_Lab/workspace/aaron/strawberry/Royal_Royce 
#minimap2 -x asm5 -t 96 farr1.fa ragtag.scaffolds.curated.reorientated.fasta > RR.GP33.minimap2.paf 

## Chromosome labelling was based on 'Camarosa'
## Renaming them figures in a more meaningful way, i.e. based on subgenomes

GP33_RR <- read_paf("RR.GP33.minimap2.paf")%>% 
  mutate(
    qname = str_replace_all( qname, pattern = c("Fvb1-1_RagTag" = "Fchil1-B1", "Fvb1-2_RagTag" = "Fchil1-Bi", "Fvb1-3_RagTag" = "Fchil1-B2", "Fvb1-4_RagTag" ="Fchil1-Av",
                                                "Fvb2-1_RagTag" = "Fchil2-B1", "Fvb2-2_RagTag" = "Fchil2-Av", "Fvb2-3_RagTag" = "Fchil2-B2", "Fvb2-4_RagTag" ="Fchil2-Bi",
                                                "Fvb3-1_RagTag" = "Fchil3-B1", "Fvb3-2_RagTag" = "Fchil3-Bi", "Fvb3-3_RagTag" = "Fchil3-B2", "Fvb3-4_RagTag" ="Fchil3-Av",
                                                "Fvb4-1_RagTag" = "Fchil4-B1", "Fvb4-2_RagTag" = "Fchil4-B2", "Fvb4-3_RagTag" = "Fchil4-Av", "Fvb4-4_RagTag" ="Fchil4-Bi",
                                                "Fvb5-1_RagTag" = "Fchil5-Av", "Fvb5-2_RagTag" = "Fchil5-B2", "Fvb5-3_RagTag" = "Fchil5-Bi", "Fvb5-4_RagTag" ="Fchil5-B1",
                                                "Fvb6-1_RagTag" = "Fchil6-Av", "Fvb6-2_RagTag" = "Fchil6-B1", "Fvb6-3_RagTag" = "Fchil6-Bi", "Fvb6-4_RagTag" ="Fchil6-B2",
                                                "Fvb7-1_RagTag" = "Fchil7-B2", "Fvb7-2_RagTag" = "Fchil7-Av", "Fvb7-3_RagTag" = "Fchil7-Bi", "Fvb7-4_RagTag" ="Fchil7-B1")))


B1 <- c("Fvb1-1", "Fvb2-1", "Fvb3-1", "Fvb4-1", "Fvb5-4", "Fvb6-2", "Fvb7-4")
Bi <- c("Fvb1-2", "Fvb2-4", "Fvb3-2", "Fvb4-4", "Fvb5-3", "Fvb6-3", "Fvb7-3")
B2 <- c("Fvb1-3", "Fvb2-3", "Fvb3-3", "Fvb4-2", "Fvb5-2", "Fvb6-4", "Fvb7-1")
Av <- c("Fvb1-4", "Fvb2-2", "Fvb3-4", "Fvb4-3", "Fvb5-1", "Fvb6-1", "Fvb7-2")


to_keep <- list(
  c("Fchil1-Av", "Fchil1-B1", "Fchil1-Bi", "Fchil1-B2",
    "Fchil2-Av", "Fchil2-B1", "Fchil2-Bi", "Fchil2-B2",
    "Fchil3-Av", "Fchil3-B1", "Fchil3-Bi", "Fchil3-B2",
    "Fchil4-Av", "Fchil4-B1", "Fchil4-Bi", "Fchil4-B2",
    "Fchil5-Av", "Fchil5-B1", "Fchil5-Bi", "Fchil5-B2",
    "Fchil6-Av", "Fchil6-B1", "Fchil6-Bi", "Fchil6-B2",
    "Fchil7-Av", "Fchil7-B1", "Fchil7-Bi", "Fchil7-B2"),
  c("chr_1A", "chr_1D", "chr_1B", "chr_1C",
    "chr_2A", "chr_2C", "chr_2B", "chr_2D",
    "chr_3A", "chr_3D", "chr_3B", "chr_3C",
    "chr_4A", "chr_4D", "chr_4B", "chr_4C",
    "chr_5A", "chr_5C", "chr_5B", "chr_5D",
    "chr_6A", "chr_6C", "chr_6B", "chr_6D",
    "chr_7A", "chr_7D", "chr_7B", "chr_7C")
)

dot_GP33_RR <- dotplot1(GP33_RR, label_seqs=TRUE, order_by="provided", ordering=to_keep,
                        xlab = "F. chiloensis", ylab = "F. x ananassa Royal Royce")  +
  theme_bw()
dot_GP33_RR
ggsave(dot_GP33_RR,filename = "/home/caroline/Documents/Caroline_postdoc_OSU/pacbio_chiloensis/2021_Phil/GP33_RR.jpeg",
       #ggsave(dot_GP33_RR,filename = "/home/caroline/Documents/Caroline_postdoc_OSU/pacbio_chiloensis/2021_Phil/GP33_RR.pdf",
       width=30)
#width=20)

#For main Fig - subset = Group 6 homeologs against vesca, a subgenome = a color
#"Av" = "red", "B1" = "#CCBB44", "B2" = "#228833", "Bi" = "#4477AA"
#cannot sirectly do a facet with pafr - need to plot separately in different colors then merge them together at the end

to_keep_6Av <- list(
  c("Fchil6-Av"),
  c("Fvb6")
)

dot_GP33_vesca_6Av <- dotplot1(ali_ragtag_curated_scaff_ves, label_seqs=FALSE, order_by="provided", ordering=to_keep_6Av,
                           xlab = "Fchil6-Av", ylab = "Fvb6", alignment_colour = "red")  +
  theme_bw() #+
  #theme(legend.position = "none", axis.text.y = element_blank(), 
  #      axis.title.y=element_blank(), axis.ticks.y = element_blank())
  #theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
dot_GP33_vesca_6Av

to_keep_6B1 <- list(
  c("Fchil6-B1"),
  c("Fvb6")
)

dot_GP33_vesca_6B1 <- dotplot1(ali_ragtag_curated_scaff_ves, label_seqs=FALSE, order_by="provided", ordering=to_keep_6B1,
                               xlab = "Fchil6-B1", ylab = "Fvb6", alignment_colour = "#CCBB44")  +
  theme_bw()+
  theme(legend.position = "none", axis.text.y = element_blank(), 
        axis.title.y=element_blank(), axis.ticks.y = element_blank())#+
  #theme(legend.position = "none", axis.text.y = element_text(colour="white"), 
  #      axis.title.y=element_text(colour="white"), axis.ticks.y =element_line(colour="white"))
        #plot.margin = unit(c(0, 0, 0, 0), "cm"))
dot_GP33_vesca_6B1

to_keep_6Bi <- list(
  c("Fchil6-Bi"),
  c("Fvb6")
)

dot_GP33_vesca_6Bi <- dotplot1(ali_ragtag_curated_scaff_ves, label_seqs=FALSE, order_by="provided", ordering=to_keep_6Bi,
                               xlab = "Fchil6-Bi", ylab = "Fvb6", alignment_colour = "#4477AA")  +
  theme_bw() +
  theme(legend.position = "none", axis.text.y = element_blank(), 
        axis.title.y=element_blank(), axis.ticks.y = element_blank())
        #plot.margin = unit(c(0, 0, 0, 0), "cm"))
dot_GP33_vesca_6Bi

to_keep_6B2 <- list(
  c("Fchil6-B2"),
  c("Fvb6")
)

dot_GP33_vesca_6B2 <- dotplot1(ali_ragtag_curated_scaff_ves, label_seqs=FALSE, order_by="provided", ordering=to_keep_6B2,
                               xlab = "Fchil6-B2", ylab = "Fvb6", alignment_colour = "#228833")  +
  theme_bw() +
  theme(legend.position = "none", axis.text.y = element_blank(), 
        axis.title.y=element_blank(), axis.ticks.y = element_blank())
       #plot.margin = unit(c(0, 0, 0, 0), "cm"))
dot_GP33_vesca_6B2

# Saving the chr6 homeologs dotplots all together
#first_line = plot_grid(dot_GP33_vesca_6Av, dot_GP33_vesca_6B1, dot_GP33_vesca_6B2, dot_GP33_vesca_6Bi, ncol = 4, align = "hv", axis = "lr") 
#first_line
#jpeg(filename = "/home/caroline/Documents/Caroline_postdoc_OSU/pacbio_chiloensis/2021_Phil/Fig_6homeologs_dotplot_GP33_vesca1.jpeg")
ggsave(
  wrap_plots(dot_GP33_vesca_6Av, dot_GP33_vesca_6B1, dot_GP33_vesca_6B2, dot_GP33_vesca_6Bi, nrow=1),
  filename = "/home/caroline/Documents/Caroline_postdoc_OSU/pacbio_chiloensis/2021_Phil/Fig_6homeologs_dotplot_GP33_vesca2.jpeg")
#dev.off()
#ggsave(first_line,filename = "/home/caroline/Documents/Caroline_postdoc_OSU/pacbio_chiloensis/2021_Phil/Fig_6homeologs_dotplot_GP33_vesca.jpeg")

