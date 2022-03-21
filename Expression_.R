#Extract the sig. MHC genes from normalized counts matrix 
HLA_A <- plotCounts(dds_DESeq, gene = "", intgroup = "MHC_I", 
                    returnData = T) 
HLA_A <- HLA_A %>%
  rownames_to_column(var = "HLA_A")
HLA_B<- plotCounts(dds_DESeq, gene = "ENSG00000234745", intgroup = "MHC_I", 
                   returnData = T) 
HLA_B <- HLA_B %>%
  rownames_to_column(var = "HLA_B")
HLA_C <- plotCounts(dds_DESeq, gene = "ENSG00000204525", intgroup = "MHC_I", 
                    returnData = T)
HLA_C <- HLA_C %>%
  rownames_to_column(var = "HLA_C")
PDCD1 <- plotCounts(dds_DESeq, gene = "ENSG00000188389", intgroup = "MHC_I", 
                    returnData = T)
PDCD1 <- PDCD1 %>%
  rownames_to_column(var = "PDCD1")
CD274 <- plotCounts(dds_DESeq, gene = "ENSG00000120217", intgroup = "MHC_I", 
                    returnData = T)
CD274 <- CD274 %>%
  rownames_to_column(var = "CD274")
CTLA4 <- plotCounts(dds_DESeq, gene = "ENSG00000163599", intgroup = "MHC_II", 
                    returnData = T)
CTLA4 <- CTLA4 %>%
  rownames_to_column(var = "CTLA4")
LAG3 <- plotCounts(dds_DESeq, gene = "ENSG00000089692", intgroup = "MHC_II", 
                   returnData = T)
LAG3 <- LAG3 %>%
  rownames_to_column(var = "LAG3")
IDO1 <- plotCounts(dds_DESeq, gene = "ENSG00000131203", intgroup = "MHC_II", 
                   returnData = T)
IDO1 <- IDO1 %>%
  rownames_to_column(var = "IDO1")
IFNG <- plotCounts(dds_DESeq, gene = "ENSG00000111537", intgroup = "MHC_II", 
                   returnData = T)
IFNG <- IFNG %>%
  rownames_to_column(var = "IFNG")

HAVCR2 <- plotCounts(dds_DESeq, gene = "ENSG00000135077", intgroup = "MHC_II", 
                     returnData = T)
HAVCR2 <- HAVCR2 %>%
  rownames_to_column(var = "HAVCR2")

HLA_DQA1 <- plotCounts(dds_DESeq, gene = "ENSG00000196735", intgroup = "MHC_II", 
                       returnData = T)
HLA_DQA1 <- HLA_DQA1 %>%
  rownames_to_column(var = "HLA_DQA1")

HLA_DQA2 <- plotCounts(dds_DESeq, gene = "ENSG00000237541", intgroup = "MHC_II", 
                       returnData = T) 
HLA_DQA2 <- HLA_DQA2 %>%
  rownames_to_column(var = "HLA_DQA2")

HLA_DQB1 <- plotCounts(dds_DESeq, gene = "ENSG00000179344", intgroup = "MHC_II", 
                       returnData = T)
HLA_DQB1 <- HLA_DQB1 %>%
  rownames_to_column(var = "HLA_DQB1")

HLA_DQB2 <- plotCounts(dds_DESeq, gene = "ENSG00000232629", intgroup = "MHC_II", 
                       returnData = T) 
HLA_DQB2 <- HLA_DQB2 %>%
  rownames_to_column(var = "HLA_DQB2")

HLA_DRA <- plotCounts(dds_DESeq, gene = "ENSG00000204287", intgroup = "MHC_II", 
                      returnData = T) 
HLA_DRA <- HLA_DRA %>%
  rownames_to_column(var = "HLA_DRA")

HLA_DRB1 <- plotCounts(dds_DESeq, gene = "ENSG00000196126", intgroup = "MHC_II", 
                       returnData = T) 
HLA_DRB1 <- HLA_DRB1 %>%
  rownames_to_column(var = "HLA_DRB1")

HLA_DRB5 <- plotCounts(dds_DESeq, gene = "ENSG00000198502", intgroup = "MHC_II", 
                       returnData = T)
HLA_DRB5 <- HLA_DRB5 %>%
  rownames_to_column(var = "HLA_DRB5")

HLA_DRB6 <- plotCounts(dds_DESeq, gene = "ENSG00000229391", intgroup = "MHC_II", 
                       returnData = T)
HLA_DRB6 <- HLA_DRB6 %>%
  rownames_to_column(var = "HLA_DRB6")

HLA_DPA1 <- plotCounts(dds_DESeq, gene = "ENSG00000231389", intgroup = "MHC_II", 
                       returnData = T)
HLA_DPA1 <- HLA_DPA1 %>%
  rownames_to_column(var = "HLA_DPA1")

HLA_DPB1 <- plotCounts(dds_DESeq, gene = "ENSG00000223865", intgroup = "MHC_II", 
                       returnData = T)
HLA_DPB1 <- HLA_DPB1 %>%
  rownames_to_column(var = "HLA_DPB1")

#######################################

tabela_exp_MHC_I_IC <- inner_join(HLA_A, HLA_B, by = c("HLA_A" = "HLA_B")) %>%
  left_join(HLA_C, by = c("HLA_A" = "HLA_C")) %>%
  left_join(PDCD1, by = c("HLA_A" = "PDCD1")) %>%
  left_join(CD274, by = c("HLA_A" = "CD274")) %>%
  left_join(CTLA4, by = c("HLA_A" = "CTLA4")) %>%
  left_join(LAG3, by = c("HLA_A" = "LAG3")) %>%
  left_join(IDO1, by = c("HLA_A" = "IDO1")) %>%
  left_join(IFNG, by = c("HLA_A" = "IFNG")) %>%
  left_join(HAVCR2, by = c("HLA_A" = "HAVCR2")) 

tabela_exp_MHC_II_IC <- inner_join(HLA_DQA1, HLA_DQA2, by = c("HLA_DQA1" = "HLA_DQA2")) %>%
  left_join(HLA_DQB1, by = c("HLA_DQA1" = "HLA_DQB1")) %>%
  left_join(HLA_DQB2, by = c("HLA_DQA1" = "HLA_DQB2")) %>%
  left_join(HLA_DRA, by = c("HLA_DQA1" = "HLA_DRA")) %>%
  left_join(HLA_DRB1, by = c("HLA_DQA1" = "HLA_DRB1")) %>%
  left_join(HLA_DRB5, by = c("HLA_DQA1" = "HLA_DRB5")) %>%
  left_join(HLA_DRB6, by = c("HLA_DQA1" = "HLA_DRB6")) %>%
  left_join(HLA_DPA1, by = c("HLA_DQA1" = "HLA_DPA1")) %>%
  left_join(HLA_DPB1, by = c("HLA_DQA1" = "HLA_DPB1")) %>%
  left_join(PDCD1, by = c("HLA_DQA1" = "PDCD1")) %>%
  left_join(CD274, by = c("HLA_DQA1" = "CD274")) %>%
  left_join(CTLA4, by = c("HLA_DQA1" = "CTLA4")) %>%
  left_join(LAG3, by = c("HLA_DQA1" = "LAG3")) %>%
  left_join(IDO1, by = c("HLA_DQA1" = "IDO1")) %>%
  left_join(IFNG, by = c("HLA_DQA1" = "IFNG")) %>%
  left_join(HAVCR2, by = c("HLA_DQA1" = "HAVCR2"))

tabela_exp_MHC_I_IC[ , c(5, 7, 9, 11, 13, 15, 17, 19, 21)] <- NULL 

tabela_exp_MHC_II_IC[ , c(5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 
                          29, 31, 33, 35)] <- NULL
names <- c("Patient_ID", "count_HLA_DQA1", "MHC_II_status", "count_HLA_DQA2",
           "count_HLA_DQB1", "count_HLA_DQB2", "count_HLA_DRA", "count_HLA_DRB1",
           "count_HLA_DRB5", "count_HLA_DRB6", "count_HLA_DPA1", "count_HLA_DPB1",
           "count_PDCD1", "count_CD274", "count_CTLA4", "count_LAG3", "count_IDO1", 
           "count_IFNG", "count_HAVCR2")

colnames(tabela_exp_MHC_II_IC) <- names
tabela_exp_MHC_II_IC <- tabela_exp_MHC_II_IC[, c("Patient_ID", "MHC_II_status", "count_HLA_DQA1", "count_HLA_DQA2",
                                               "count_HLA_DQB1", "count_HLA_DQB2", "count_HLA_DRA", "count_HLA_DRB1",
                                               "count_HLA_DRB5", "count_HLA_DRB6", "count_HLA_DPA1", "count_HLA_DPB1",
                                               "count_PDCD1", "count_CD274", "count_CTLA4", "count_LAG3", "count_IDO1", 
                                               "count_IFNG", "count_HAVCR2")]

tabela_exp_cor <- round(cor(tabela_exp_MHC_II_IC[ , 3:19]), 4)
tabela_exp_cor_pvalue <- cor_pmat(tabela_exp_MHC_II_IC[ , 3:19])
pheatmap(tabela_exp_cor)

write.table(tabela_exp_cor, "mCRPC/Results/MHC/MHC_II/table_correlation_MHC_II.txt",sep = "\t", 
            row.names = T, col.names = T)

library(ggpubr)
ggplot(tabela_exp, mapping = aes(x = count_CDK12, y = count_JAK1 , color = CDK12_status)) +
  geom_point() +
  geom_smooth(method = lm, se = F) +
  theme_bw() +
  stat_cor(method = "pearson", label.x = 3) +
  scale_y_log10() +
  scale_x_log10() +
  ggtitle("Pearson Correlation: CDK12 vs JAK1")

library(ggcorrplot)
install.packages('ggcorrplot')

ggcorrplot(tabela_exp_cor, hc.order = TRUE, outline.col = "white", p.mat = tabela_exp_cor_pvalue, 
           title = "Pearson Correlation", lab = F, pch = T, pch.cex = T) 


ggcorrplot(tabela_exp_cor, hc.order = TRUE, type = "lower",
           outline.col = "white")
