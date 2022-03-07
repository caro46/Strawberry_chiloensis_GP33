# Creating Figure for the SDR genes phylogeny
# This was used to produce the main figure 6 of Cauret et al. 2022
# Chromosome-scale assembly with a phased sex-determining region resolves features of early Z and W chromosome differentiation in a wild octoploid strawberry

# Loading libraries
library(ggtree)
library(ape)
library(treeio)
library(dplyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(ggrepel)

# Setting working directoy
setwd("/home/caroline/Documents/Caroline_postdoc_OSU/pacbio_chiloensis/SDR_phylo")

# Loading some custom themes for plotting
source("~/Documents/Caroline_postdoc_OSU/github_repo/Strawberry_chiloensis_GP33/R_plots_themes.R")

# Reading the tree files (outputs from iqtree)

## Below were attempts with the files AL saved from figtree - did not work properly
#treeRPPO <-read.tree("SDR_RPP0.mafft.fsa.treefile")
#treeRPPO <-read.nexus("SDR_RPP0W_mafft.fasta.treefile")
#treeRPPO <-read.beast("SDR_RPP0W_mafft.fasta.treefile")
#treeGLU <-read.tree("SDR_glucan_endo-1.mafft.fsa.treefile.min")
#treeGLU <-read.nexus("SDR_glucan_endo-1.mafft.fsa.treefile")
#ggtree(treeGLU)
#treeGME <-read.tree("SDR_GME_complete.mafft.fsa.treefile.min")
#treeGME <-read.nexus("SDR_GME_complete.mafft.fsa.treefile")

## Original iqtree output (did not go through figtree first)
treeRPPO <-read.tree("SDR_RPPO_iqtree_out.tre")
treeGLU <- read.tree("SDR_glucan_iqtree_out.tre")
treeGME <-read.tree("SDR_GME_iqtree_out.tre")

# Obtaining dataframe from tree
treeRPPO_df <-as_tibble(treeRPPO)
treeGLU_df <-as_tibble(treeGLU)
treeGME_df <-as_tibble(treeGME)
#View(treeRPPO_df)

# Adding columns to be used for coloring and shaping and keeping only useful columns
treeRPPO_df <-treeRPPO_df %>% 
  mutate(species = case_when( str_detect(label, regex("FxaC")) ~ "camarosa",
                              str_detect(label, regex("Fchil")) ~ "chiloensis",
                              str_detect(label, regex("FvH4")) ~ "vesca",
                              str_detect(label, regex("SDR")) ~ "chiloensis",
                              TRUE ~ "unlabeled")
         ) %>%
  mutate(subgenome = case_when( str_detect(label, regex("25g44180")) ~ "Av",
                              str_detect(label, regex("7g24570")) ~ "Av",
                              str_detect(label, regex("26g13040")) ~ "Bi",
                              str_detect(label, regex("27G_18678957")) ~ "B2",
                              str_detect(label, regex("SDR")) ~ "B2",
                              str_detect(label, regex("28g09761")) ~ "B1",
                              #adding subg info to "sub-clade"
                              str_detect(label, regex("25g43640")) ~ "Av",
                              str_detect(label, regex("7g24170")) ~ "Av",
                              str_detect(label, regex("26g13660")) ~ "Bi",
                              str_detect(label, regex("27g39570")) ~ "B2",
                              str_detect(label, regex("28g10250")) ~ "B1",
                              TRUE ~ "other chromosome")
  ) %>%
  select(label, species, subgenome)
# Doing a full join between the tree and othe rdata is in theory an option but froze the computer each time
# using %<+% below works perfectly
#tree_full_info_RPPO <-full_join(treeRPPO, tree_df_sub, by="label")
#tree_full_info_RPPO <-full_join(treeRPPO, treeRPPO_df, by="label")

#png("RPPO_phylo_try.png")
#cols <- c("Av" = "red", "Bi" = "#4477AA", "B1" = "#CCBB44", "B2" = "#228833", "other chromosome" = "black")
cols <- c("Av" = "red", "Bi" = "#4477AA", "B1" = "#CCBB44", "B2" = "#228833")
shapes_sp <- c("camarosa" = 18, "chiloensis" = 17, "vesca" = 15)

# First look at the unrooted tree
ggtree(treeRPPO) %<+% treeRPPO_df + geom_tippoint(aes(shape=species, color = subgenome), size = 3) +
  geom_tiplab() + 
  #geom_text(aes(x=branch, label=label), 
  #          vjust=-.03, size=1.8) +
  #geom_nodelab(aes(label = node))+
  #geom_text2(aes(subset=!isTip, label=label), hjust=-.2, size=4)+
  scale_color_manual(values = cols) +
  #geom_text2(aes(subset=!isTip, label=node))+
  geom_text2(aes(label=node))+
  #scale_fill_manual(values = cols) +
  scale_shape_manual(values = shapes_sp)

# Rooting the tree
RPPO_tree_root <- ape::root(treeRPPO, node = 62, edgelabel = TRUE)
#pdf("RPPO_phylo_try_root.pdf")
#jpeg("RPPO_phylo_try_root.jpeg")
RPPO_root_tree <- ggtree(RPPO_tree_root) %<+% treeRPPO_df + geom_tippoint(aes(shape=species, color = subgenome), size = 3) +
  #geom_tiplab() + 
  #geom_text2(aes(label=node)) +
  #geom_nodelab()+
  scale_color_manual(values = cols) +
  scale_shape_manual(values = shapes_sp) +
  geom_cladelabel(node=2, label="SDR", geom='label', 
                  color='black', fontsize=3, fill='#228833', offset = -0.003) +
  geom_cladelabel(node=1, label="Fchil7-B2",color='black', fontsize=3, offset = -0.003) +
  #geom_point2(size = 3, shape = 21, fill= "black", aes(subset=!isTip & !is.na(as.numeric(label)) & as.numeric(label) >= 80)) +
  geom_point2(size = 3, shape = 21, fill= "white", aes(subset=!isTip & !is.na(as.numeric(label)) & as.numeric(label) < 80)) +
  #theme_tree2() +
  #xlim(0, 0.2) + #needed for subset 
  #geom_treescale(x=0.14, y=30)+
  theme(legend.position = "none")
  #theme(legend.position='top', 
  #      legend.justification='left',
  #      legend.direction='vertical')
#dev.off()
#RPPO_root_tree %>% collapse(node = 40, 'min')%>% 
#  collapse(node = 40, 'min')

# Subsetting and adding a scale
RPPO_root_tree_sub <- RPPO_root_tree + xlim(0, 0.21) + geom_treescale(x=0.14, y=30) +#geom_treescale(x=0.18, y=18)+  
  geom_hilight(node = 40, fill = "grey", alpha = 0.2) +
  annotate("text", label = "duplicate", x=0.18, y=34, size = 3)
SDR_clade_RPP <- viewClade(RPPO_root_tree_sub,50)#%>% collapse(node = 39, 'min')

#dev.off()
#p <- p + geom_text(aes(label = bootstrap), hjust = 1, vjust = -0.4, size = 3) + geom_nodelab(aes(label = nodesupport)) # specify your node label here, looks like BP
#View(treeRPPO_df)
#ggtree(treeRPPO) + 
#  geom_tippoint(data = treeRPPO_df, aes(color=subgenome), size=1.5)
#ggtree(tree_full_info_RPPO) + #freezes whole computer
#  geom_tippoint(aes(color=subgenome), size=1.5)
#dev.off()

# Same for GLU
# Adding columns to be used for coloring and shaping and keeping only useful columns
treeGLU_df <-treeGLU_df %>%
  mutate(species = case_when( str_detect(label, regex("FxaC")) ~ "camarosa",
                              str_detect(label, regex("Fchil")) ~ "chiloensis",
                              str_detect(label, regex("FvH4")) ~ "vesca",
                              str_detect(label, regex("SDR")) ~ "chiloensis",
                              TRUE ~ "unlabeled")
  ) %>%
  mutate(subgenome = case_when( str_detect(label, regex("21g36850")) ~ "Av",
                                str_detect(label, regex("6g24680")) ~ "Av",
                                str_detect(label, regex("5g06210")) ~ "Av",
                                str_detect(label, regex("18g39760")) ~ "Bi",
                                str_detect(label, regex("24g31020")) ~ "B2",
                                str_detect(label, regex("SDR")) ~ "B1",
                                str_detect(label, regex("23g15920")) ~ "B1",
                                TRUE ~ "other chromosome")
  ) %>%
  select(label, species, subgenome)

#View(treeGLU_df)

# First look at the unrooted tree
ggtree(treeGLU) %<+% treeGLU_df + geom_tippoint(aes(shape=species, color = subgenome), size = 3) +
  geom_tiplab() + 
  geom_text2(aes(label=node)) +
  #geom_nodelab(size = 7, color= "red")+
  #geom_nodelab(geom='label',size = 7, color= "red", aes(label=label, subset=!is.na(as.numeric(label)) & as.numeric(label) < 80))+
  scale_color_manual(values = cols) 

# Rooting
GLU_tree_root <- ape::root(treeGLU, node = 28, edgelabel = TRUE)
GLU_tree_root <- ggtree(GLU_tree_root) %<+% treeGLU_df + geom_tippoint(aes(shape=species, color = subgenome), size = 3) +
  #geom_tiplab() + 
  #geom_text2(aes(label=node)) +
  scale_color_manual(values = cols) +
  scale_shape_manual(values = shapes_sp) +
  #geom_cladelabel(node=1, label="SDR", geom='label', 
  #                color='black', fontsize=3, fill='#CCBB44') +
  geom_cladelabel(node=1, label="SDR", geom='label', 
                  color='black', fontsize=3, fill='#CCBB44', offset = -0.009) +
  geom_cladelabel(node=21, label="Fchil6-B1",color='black', fontsize=3, offset = -0.009) +
  #location - not needed anymore
  # geom_cladelabel(node=18, label="Fvb7", 
  #                 color='black', fontsize=2) +
  # geom_cladelabel(node=17, label="Fvb7-2", 
  #                 color='black', fontsize=2) +
  # geom_cladelabel(node=20, label="Fchil-B2", 
  #                 color='black', fontsize=2) +
  # geom_cladelabel(node=19, label="Fvb7-1", 
  #                 color='black', fontsize=2) +
  # geom_cladelabel(node=21, label="Fvb7-4", 
  #                 color='black', fontsize=2) +
  #geom_point2(size = 3, shape = 21, fill= "black", aes(subset=!isTip & !is.na(as.numeric(label)) & as.numeric(label) >= 80)) +
  geom_point2(size = 3, shape = 21, fill= "white", aes(subset=!isTip & !is.na(as.numeric(label)) & as.numeric(label) < 80)) +
  #theme_tree2() +
  theme(legend.position = "none") #+
  #xlim(0, 0.4) +
  #geom_nodelab(geom='label',size = 3, color= "black", aes(label=label, subset=!is.na(as.numeric(label)) & as.numeric(label) < 80),alpha=0.5,
  #             label.size = 0)
  #geom_label_repel(aes(label=label, subset=!is.na(as.numeric(label)) & as.numeric(label) < 80))

GLU_tree_root
GLU_root_tree_sub <- GLU_tree_root + geom_treescale(x=0.335, y=12.5) + xlim(0.31, 0.35)  + ylim(12,21)+#xlim(0.3, 0.4)  + ylim(12,21) +#geom_treescale(x=0.28, y=20)+#geom_treescale(x=0.31, y=20)
#geom_cladelab(node=27, label="Group 5 homeologs", align=TRUE, barsize=1.5)
  geom_hilight(node = 27, fill = "grey", alpha = 0.2) +
  annotate("text", label = "chromosome 5 homeologs", x=0.326, y=15.25, size = 3)+
  annotate("segment", x = 0.314, y = 15.88, 
           #xend = 0.317, yend = 15.88, color = "black") 
           xend = 0.3166, yend = 15.88, color = "black") 

  #geom_hilight(node = 40, fill = "grey", alpha = 0.2)
#SDR_clade_GLU <- viewClade(GLU_tree_root,26)  #not working well because long branch

# tree_full_info_GLU <-full_join(treeGLU, treeGLU_df, by="label")
# ggtree(tree_full_info_GLU)

# Same for GME
# Adding columns to be used for coloring and shaping and keeping only useful columns
treeGME_df <-treeGME_df %>%
  mutate(species = case_when( str_detect(label, regex("FxaC")) ~ "camarosa",
                              str_detect(label, regex("Fchil")) ~ "chiloensis",
                              str_detect(label, regex("SDR")) ~ "chiloensis",
                              str_detect(label, regex("FvH4")) ~ "vesca",
                              TRUE ~ "unlabeled")
  ) %>%
  mutate(subgenome = case_when( str_detect(label, regex("21g67460")) ~ "Av",
                                str_detect(label, regex("6g02980")) ~ "Av",
                                str_detect(label, regex("22g66610")) ~ "Bi",
                                str_detect(label, regex("22g67740")) ~ "Bi",
                                str_detect(label, regex("24g07620")) ~ "B2",
                                str_detect(label, regex("SDR")) ~ "B2",
                                str_detect(label, regex("23g41100")) ~ "B1",
                                TRUE ~ "other chromosome")
  )%>%
  select(label, species, subgenome)

#View(treeGME_df)

# First look at the unrooted tree
GME_tree_unroot <- ggtree(treeGME) %<+% treeGME_df + geom_tippoint(aes(shape=species, color = subgenome), size = 3) +
  geom_tiplab() + 
  geom_text2(aes(label=node)) +
  scale_color_manual(values = cols) 

# Rooting
GME_tree_root <- ape::root(treeGME, node = 19, edgelabel = TRUE)
GME_tree_root <- ggtree(GME_tree_root) %<+% treeGME_df + geom_tippoint(aes(shape=species, color = subgenome), size = 3) +
  #geom_tiplab() +
  #geom_nodelab(geom='label',size = 3, color= "black", aes(label=label, subset=!is.na(as.numeric(label))),alpha=0.5,
  #             label.size = 0)+
  #geom_text2(aes(label=node)) +
  scale_color_manual(values = cols, guide = "textcolourguide") +
  scale_shape_manual(values = shapes_sp, labels = c("Camarosa", "F. chiloensis", "F. vesca"), guide =
                       guide_legend(label.theme = element_text(angle = 0, face = "italic", size = 8))) +
  geom_cladelabel(node=1, label="SDR", geom='label', 
                  color='black', fontsize=3, fill='#228833') +
  geom_cladelabel(node=11, label="Fchil6-B2",color='black', fontsize=3, offset = 0.0005) +
  geom_cladelabel(node=10, label="Fvb6-4",color='black', fontsize=3, offset = 0.0005) +
  #theme_tree2() +
  geom_treescale(x=0.005, y = 1) +
  xlim(0, 0.02) +
  #geom_point2(size = 3, shape = 21, fill= "black", aes(subset=!isTip & !is.na(as.numeric(label)) & as.numeric(label) >= 80)) +
  geom_point2(size = 3, shape = 21, fill= "white", aes(subset=!isTip & !is.na(as.numeric(label)) & as.numeric(label) < 80)) #+
  #theme(legend.position = "none") 
GME_tree_root

# Getting the legend as a separate plot 
legend_phylo <- cowplot::get_legend(GME_tree_root) 
GME_tree_root <- GME_tree_root +  theme(legend.position = "none") 

# Combining the 3 phylogenies
#SDR_genes_phylo <- plot_grid(RPPO_root_tree, GME_tree_root, GLU_tree_root, rel_heights = c(1, 1,1), ncol = 3, align = "v", labels=c('A) RPP0', ' B) GME', 'C) glucan'))
#ggsave(SDR_genes_phylo,filename = "SDR_genes_phylo.jpeg")
#ggsave(SDR_genes_phylo,filename = "SDR_genes_phylo.pdf", width = 14 , height = 7, units = "in")

#SDR_genes_phylo <- plot_grid(RPPO_root_tree, GME_tree_root, GLU_tree_root, legend_phylo,rel_widths = c(3, 3,3,1), ncol = 4, align = "v", labels=c('A) RPP0', ' B) GME', 'C) glucan'))
#SDR_genes_phylo
#ggsave(SDR_genes_phylo,filename = "SDR_genes_phylo_general_legend.pdf", width = 16 , height = 6, units = "in")
#ggsave(SDR_genes_phylo,filename = "SDR_genes_phylo_general_legend_support.pdf", width = 16 , height = 6, units = "in")

# Best results
SDR_genes_phylo_sub <- plot_grid(SDR_clade_RPP, GME_tree_root, GLU_root_tree_sub, legend_phylo,rel_widths = c(3, 3,3,1), ncol = 4, align = "v", labels=c('A) RPP0', ' B) GME', 'C) glucan'))
SDR_genes_phylo_sub
#ggsave(SDR_genes_phylo_sub,filename = "SDR_genes_phylo_general_legend_support_subset.pdf", width = 16 , height = 6, units = "in")
save_plot("SDR_genes_phylo_general_legend_support_subset.pdf", SDR_genes_phylo_sub, ncol = 4, base_asp = 1.1)
#save_plot("file.jpeg", SDR_genes_phylo_sub, ncol = 4, base_asp = 1.1)


