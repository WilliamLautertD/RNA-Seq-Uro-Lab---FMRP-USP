library(ggplot2)
library(ggsignif)
library(ggpubr)
library(reshape2)

gep <- read.csv(' ', 
                header = T, sep = ',')
  
metadata <- read.csv('', header = T,
                      sep = "\t")
rownames(metadata) <- metadata$Mixture

CDK12_GENOME_mCRPC <- merge(gep, metadata, by =  "Mixture")

CDK12_GENOME_mCRPC$MHC_I <- factor(CDK12_GENOME_mCRPC$MHC_I, 
                                          levels = c("CDK12_Hom", "CDK12_Het", "CDK12_WT"))
CDK12_GENOME_mCRPC$MHC_II <- factor(CDK12_GENOME_mCRPC$MHC_II, 
                                   levels = c("High", "Low"))

instab_frame <- CDK12_GENOME_mCRPC[,c( "",
                                       "B.cells.naive",
                                       "B.cells.memory",
                                       "Plasma.cells",
                                       "T.cells.CD8",
                                       "T.cells.CD4.naive",
                                       "T.cells.CD4.memory.resting",
                                       "T.cells.CD4.memory.activated",
                                       "T.cells.follicular.helper",
                                       "T.cells.regulatory..Tregs.",
                                       "T.cells.gamma.delta",
                                       "NK.cells.resting",
                                       "NK.cells.activated",
                                       "Monocytes",
                                       "Macrophages.M0",
                                       "Macrophages.M1",
                                       "Macrophages.M2",
                                       "Dendritic.cells.resting",
                                       "Dendritic.cells.activated",
                                       "Mast.cells.resting",
                                       "Mast.cells.activated",
                                       "Eosinophils",
                                       "Neutrophils")]



instabstat <- melt(instab_frame, id.var = '')
ggplot(data = na.omit(instabstat), aes(x = variable, y = value)) +
  ggtitle("") +
  scale_fill_manual(values = c("red2", "dodgerblue2", "forestgreen", "lightblue")) +
  theme_minimal() +
  geom_point(position = position_jitter(w = 0.3,h = 0), alpha=0.1) +
  theme(axis.text = element_text(size = 10, face = 'bold'),
        axis.title = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        plot.title = element_text(hjust = 0, size = 10),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13))+
  geom_boxplot(aes(fill = ), outlier.size = 0.07, alpha = 0.7, size = 0.5)+
  scale_x_discrete(name = "")+
  scale_y_continuous(name = "CIBERSORT relative abundance")+
  stat_compare_means(aes(group = ), paired = F, method = "t.test", label = "p.signif", hide.ns = T)
