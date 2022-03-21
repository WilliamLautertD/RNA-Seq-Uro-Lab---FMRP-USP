##Obtain logical vector where TRUE values denote padj values < 0.05 and 
HLA_A <- plotCounts(dds_DESeq, gene = "ENSG00000206503", intgroup = "MHC_I", 
                    returnData = T) 
HLA_B<- plotCounts(dds_DESeq, gene = "ENSG00000234745", intgroup = "MHC_I", 
                   returnData = T) 
HLA_C <- plotCounts(dds_DESeq, gene = "ENSG00000204525", intgroup = "MHC_I", 
                    returnData = T) 
HLA_DQA1 <- plotCounts(dds_DESeq, gene = "ENSG00000196735", intgroup = "MHC_II", 
                       returnData = T) 
HLA_DQA2 <- plotCounts(dds_DESeq, gene = "ENSG00000237541", intgroup = "MHC_II", 
                       returnData = T) 
HLA_DQB1 <- plotCounts(dds_DESeq, gene = "ENSG00000179344", intgroup = "MHC_II", 
                       returnData = T) 
HLA_DQB2 <- plotCounts(dds_DESeq, gene = "ENSG00000232629", intgroup = "MHC_II", 
                       returnData = T) 
HLA_DRA <- plotCounts(dds_DESeq, gene = "ENSG00000204287", intgroup = "MHC_II", 
                      returnData = T) 
HLA_DRB1 <- plotCounts(dds_DESeq, gene = "ENSG00000196126", intgroup = "MHC_II", 
                       returnData = T) 
HLA_DRB5 <- plotCounts(dds_DESeq, gene = "ENSG00000198502", intgroup = "MHC_II", 
                       returnData = T)
HLA_DRB6 <- plotCounts(dds_DESeq, gene = "ENSG00000229391", intgroup = "MHC_II", 
                       returnData = T)
HLA_DPA1 <- plotCounts(dds_DESeq, gene = "ENSG00000231389", intgroup = "MHC_II", 
                       returnData = T)
HLA_DPB1 <- plotCounts(dds_DESeq, gene = "ENSG00000223865", intgroup = "MHC_II", 
                       returnData = T)
HLA_DMA <- plotCounts(dds_DESeq, gene = "ENSG00000204257", intgroup = "MHC_II", 
                      returnData = T)
HLA_DMB <- plotCounts(dds_DESeq, gene = "ENSG00000242574", intgroup = "MHC_II", 
                      returnData = T)
summary(HLA_A)
HLA_A <- HLA_A %>% 
  mutate(threshold = count < 6456) %>%
  rownames_to_column(var = "HLA_A")

summary(HLA_B)
HLA_B <- HLA_B %>%
  mutate(threshold = count < 10671) %>%
  rownames_to_column(var = "HLA_B")

summary(HLA_C)
HLA_C <- HLA_C %>%
  mutate(threshold = count < 17305) %>%
  rownames_to_column(var = "HLA_C")

tabela_exp_MHC_I <- inner_join(HLA_A, HLA_B, by = c("HLA_A" = "HLA_B"),
                               suffix = c("HLA_A", "HLA_B")) %>%
  left_join(HLA_C, by = c("HLA_A" = "HLA_C"))

tabela_exp_MHC_I[ , c(3, 6, 9)] <- NULL
names <- c("Patient_ID", "count_HLA_A", "threshold_HLA_A", "count_HLA_B",
           "threshold_HLA_B", "count_HLA_C", "threshold_HLA_C")
colnames(tabela_exp_MHC_I) <- names

tabela_exp_MHC_I <- tabela_exp_MHC_I %>%
  mutate(group = threshold_HLA_A & threshold_HLA_B & threshold_HLA_C)
summary(tabela_exp_MHC_I)
write_csv(tabela_exp_MHC_I, "mCRPC/Results/MHC/expression/tabela_exp_MHC_I.txt", col_names = T)

summary(HLA_DQA1)
HLA_DQA1 <- HLA_DQA1 %>%
  mutate(threshold = count < 837.6) %>%
  rownames_to_column(var = "HLA_DQA1")
summary(HLA_DQA2)
HLA_DQA2 <- HLA_DQA2 %>%
  mutate(threshold = count < 52.947) %>%
  rownames_to_column(var = "HLA_DQA2")
summary(HLA_DQB1)
HLA_DQB1 <- HLA_DQB1 %>%
  mutate(threshold = count < 767.5) %>%
  rownames_to_column(var = "HLA_DQB1")
summary(HLA_DQB2)
HLA_DQB2 <- HLA_DQB2 %>% 
  mutate(threshold = count < 40.49) %>%
  rownames_to_column(var = "HLA_DQB2")
summary(HLA_DRA)
HLA_DRA <- HLA_DRA %>%
  mutate(threshold = count < 7382.5) %>%
  rownames_to_column(var = "HLA_DRA")
summary(HLA_DRB1)
HLA_DRB1 <- HLA_DRB1 %>%
  mutate(threshold = count < 2935.0) %>%
  rownames_to_column(var = "HLA_DRB1")
summary(HLA_DRB5)
HLA_DRB5 <- HLA_DRB5 %>%
  mutate(threshold = count < 587.09) %>%
  rownames_to_column(var = "HLA_DRB5")
summary(HLA_DRB6)
HLA_DRB6 <- HLA_DRB6 %>%
  mutate(threshold = count < 6.471) %>%
  rownames_to_column(var = "HLA_DRB6")
summary(HLA_DPA1)
HLA_DPA1 <- HLA_DPA1 %>%
  mutate(threshold = count < 1040.6) %>%
  rownames_to_column(var = "HLA_DPA1")
summary(HLA_DPB1)
HLA_DPB1 <- HLA_DPB1 %>%
  mutate(threshold = count < 2789.0) %>%
  rownames_to_column(var = "HLA_DPB1")

tabela_exp_MHC_II <- left_join(HLA_DQA1, HLA_DQA2, by = c("HLA_DQA1" = "HLA_DQA2"),
                               suffix = c("HLA_DQA1", "HLA_DQA2")) %>%
  left_join(HLA_DQB1, by = c("HLA_DQA1" = "HLA_DQB1")) %>%
  left_join(HLA_DQB2, by = c("HLA_DQA1" = "HLA_DQB2")) %>%
  left_join(HLA_DRA, by = c("HLA_DQA1" = "HLA_DRA")) %>%
  left_join(HLA_DRB1, by = c("HLA_DQA1" = "HLA_DRB1")) %>%
  left_join(HLA_DRB5, by = c("HLA_DQA1" = "HLA_DRB5")) %>%
  left_join(HLA_DRB6, by = c("HLA_DQA1" = "HLA_DRB6")) %>%
  left_join(HLA_DPA1, by = c("HLA_DQA1" = "HLA_DPA1")) %>%
  left_join(HLA_DPB1, by = c("HLA_DQA1" = "HLA_DPB1"))

tabela_exp_MHC_II[ , c(3, 6, 9, 12, 15, 18, 21, 24, 27, 30)] <- NULL
names <- c("Patient_ID", "count_HLA_DQA1", "threshold_HLA_DQA1", "count_HLA_DQA2",
           "threshold_HLA_DQA2", "count_HLA_DQB1", "threshold_HLA_DQB1", "count_HLA_DQB2",
           "threshold_HLA_DQB2", "count_HLA_DRA", "threshold_HLA_DRA", "count_HLA_DRB1",
           "threshold_HLA_DRB1", "count_HLA_DRB5", "threshold_HLA_DRB5", "count_HLA_DRB6", 
           "threshold_HLA_DRB6", "count_HLA_DPA1", "threshold_HLA_DPA1", "count_HLA_DPB1",
           "threshold_HLA_DPB1")
colnames(tabela_exp_MHC_II) <- names
tabela_exp_MHC_II <- tabela_exp_MHC_II %>%
  mutate(group = threshold_HLA_DQA1 & threshold_HLA_DQA2 & threshold_HLA_DQB1 & threshold_HLA_DQB2 &
           threshold_HLA_DRA & threshold_HLA_DRB1 & threshold_HLA_DRB5 & threshold_HLA_DRB6 & threshold_HLA_DPA1 &
           threshold_HLA_DPB1)
summary(tabela_exp_MHC_II)

write_csv(tabela_exp_MHC_II, "mCRPC/Results/MHC/expression/tabela_exp_MHC_II.txt", col_names = T)
