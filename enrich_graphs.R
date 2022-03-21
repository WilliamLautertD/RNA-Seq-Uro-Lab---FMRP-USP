###########
#ORA graphs 
###########
###########
#EnrichR###

e_GO_BP_Het_UP <- read_csv("enrich_GO_BP_PRAD_MHC_II_UP.txt")
e_GO_BP_Hom_UP <- read_csv("enrich_GO_BP_mCRPC_MHC_II_UP.txt")
e_GO_BP_Het_DOWN <- read_csv("enrichR/enrich_GO_BP_CDK12_Het_WT_DOWN.txt")
e_GO_BP_Hom_DOWN <- read_csv("enrichR/enrich_GO_BP_CDK12_Hom_WT_DOWN.txt")

e_GO_MF_Het_UP <- read_csv("enrichR/enrich_GO_MF_CDK12_Het_WT_UP.txt")
e_GO_MF_Hom_UP <- read_csv("enrichR/enrich_GO_MF_CDK12_Hom_WT_UP.txt")
e_GO_MF_Het_DOWN <- read_csv("enrichR/enrich_GO_MF_CDK12_Het_WT_DOWN.txt")
e_GO_MF_Hom_DOWN <- read_csv("enrichR/enrich_GO_MF_CDK12_Hom_WT_DOWN.txt")

#clusterP##
c_GO_BP_Het_UP <- read.csv("c_enrich_GO_BP_CDK12_MHC_I_High_vs_Low_UP.txt",
                           sep = "\t")
c_GO_BP_Hom_UP <- read.csv("clusterProfiler/c_enrich_GO_BP_CDK12_Hom_WT_UP.txt",
                           sep = "\t")
c_GO_BP_Het_DOWN <- read.csv("clusterProfiler/c_enrich_GO_BP_CDK12_Het_WT_DOWN.txt",
                             sep = "\t")
c_GO_BP_Hom_DOWN <- read.csv("clusterProfiler/c_enrich_GO_BP_CDK12_Hom_WT_DOWN.txt",
                             sep = "\t")

c_GO_MF_Het_UP <- read.csv("clusterProfiler/c_enrich_GO_MF_CDK12_Het_WT_UP.txt",
                           sep = "\t")
c_GO_MF_Hom_UP <- read.csv("clusterProfiler/c_enrich_GO_MF_CDK12_Hom_WT_UP.txt",
                           sep = "\t")
c_GO_MF_Het_DOWN <- read.csv("clusterProfiler/c_enrich_GO_MF_CDK12_Het_WT_DOWN.txt",
                           sep = "\t")
c_GO_MF_Hom_DOWN <- read.csv("clusterProfiler/c_enrich_GO_MF_CDK12_Hom_WT_DOWN.txt",
                           sep = "\t")
  
#####################################  
library(grid)
library(tidyverse)
library(shadowtext) 
c_GO_BP_Het_UP <- c_GO_BP_Het_UP %>%
  mutate(qscore = -log(`pvalue`, base=10)) %>%
  dplyr::arrange(qscore)

RColorBrewer::display.brewer.all()
color_up <- RColorBrewer::brewer.pal(5, name = "YlOrRd")
color_down <-  RColorBrewer::brewer.pal(5, name = "Blues")

ggplot(c_GO_BP_Het_UP[67:87,],
       aes(`qscore`, fct_reorder(Description, `qscore`), fill = `qvalue`)) +
  geom_col() +
  theme_gray() +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major.x = element_blank(),
        axis.text.y = element_blank(), 
        axis.title.y = element_blank(),
        plot.title = element_text(size = 20),
        axis.text.x = element_text(colour = 'black', size = 10),
        legend.position = c(0.9, 0.2)) +
  geom_text(mapping = aes(0, y = Description, label = Description),
            hjust = 0, colour = "black", size = 5) +
  scale_x_continuous(expand = c(0.01,0))+
  scale_y_discrete() +
  scale_fill_gradientn(colours = color_up,
                       guide = guide_colorbar(reverse = F)) +
  labs(title = "GO Molecular Function (mCRPC:CDK12-MHC_I_up)") +
  xlab("-log10(padj)")

GO_BP_TCGA_PRAD_CDK12_Hom_UP
limits = c(0, 0.03), expand = c(0.01,0)
