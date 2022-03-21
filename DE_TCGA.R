#----------------------------------------------------------------------------------------
#-------------Load Library---------------------------------------------------------------
i <- T
while(i) {
  library(DESeq2)
  library(dplyr)
  library(tidyverse)
  library(ggplot2)
  library(pheatmap)
  library(ggrepel)
  library(ggsignif)
  library(ggpubr)
  library(annotables)
  library(marray)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ensembldb)
  library(msigdbr)
  library(GSEABase)
  library(ggnewscale)
  library(vidger)
  i <- F
}
#----------------------------------------------------------------------------------------
#------------Data input (Raw counts and metadata)----------------------------------------
RNA_raw_count <- read.csv("", 
                          header = T, sep = "\t")
metadata <- read.csv("", 
                     header = T, sep = ";")

#----------------------------------------------------------------------------------------
#------------Matching both colnames and rownames of Raw_counts and Metadata tables-------
i <- T 
while(i) {
  ###---naming 'metadata' rownames---
  rownames(metadata) <- metadata$Mixture
  ###---arraging by the factor of interest---
  metadata <- metadata %>%
    arrange()
  ###---Creating the index and matching the tables
  idx <- match(rownames(metadata), colnames(RNA_raw_count))
  RNA_raw_count <- RNA_raw_count[ ,idx]
  all(rownames(metadata) == colnames(RNA_raw_count))
  ###---Return---
  if (all(rownames(metadata) == colnames(RNA_raw_count))) {
    i <- F
    print("Row names of 'Metadata' equal to Column names of 'RNA_raw_count'")
  } else {
    print("Row names of 'Metadata' is not equal to Column names of 'RNA_raw_count'")
    i <- F
  }
}

#----------------------------------------------------------------------------------------
#------------Creating the DESeq object and running DESeq2 analisys-----------------------
i <- T
while (i) {
  ###---DESeq object(Count matrix, metadata, and design formula)---
  dds_rna <- DESeqDataSetFromMatrix(RNA_raw_count, 
                                    colData = metadata, 
                                    design  = ~ )
  ###---Relevel of the design factor---
  #use the factor of interest after '$'
  dds_rna$ <- relevel(dds_rna$MHC_II, ref = "")
  
  ###---Running DESeq2---
  dds_DESeq <- DESeq(dds_rna, parallel = T)
  i <- F
}

#----------------------------------------------------------------------------------------
#------------Quality Control(QC)--------------------------------------------------------- 
###---Estimation of gene-wise dispersion and mean-variance relationship---
i <- T
while (i) {
  ###---mean---
  mean_counts <- apply(RNA_raw_count[, 1:488],1, mean) #The range in the column slot... 
                                                      #...depends of the number of samples
  ###---variantion---
  variance_counts <- apply(RNA_raw_count[, 1:488], 1, var)      
  df_int <- data.frame(mean_counts, variance_counts)
  
  ###---plot of mean-variance relationship---
  ggplot(df_int) +
    geom_point(aes(x = mean_counts, y = variance_counts)) +
    scale_x_log10(limits = c(1,1e9)) +
    scale_y_log10(limits = c(1,1e9)) +
    xlab("Mean counts per gene")     +
    ylab("Variance per gene")        +
    geom_abline(intercept = 0, slope = 1, color = "red")  +
    ggtitle("mean-variance relationship_TCGA") 
  i <- F
}

###---Plot dispersion relatioship---
plotDispEsts(dds_DESeq)

#----------------------------------------------------------------------------------------
#------------Sample-level QC:PCA and hierarchical clustering-------------------------------
i <- T
while (i) {
  vsd <- vst(dds_DESeq, blind = F)
  ###--converts the data in a DESeqTransform object to a simple 2-dimensional 
  ##data structure
  vsd_mat <- assay(vsd)
  ##Compute pairwise correlation
  vsd_cor <- cor(vsd_mat)
  i <- F
}

meta <- dplyr::select(metadata, c(MHC_II))
pheatmap(vsd_cor, annotation = meta, 
         show_rownames = F,
         show_colnames = F,
         main = "Hierarchical clustering - Correlation")

###---PCA funtion---
plotPCA(vsd, intgroup =  c("MHC_II"))

###---PCA object to use in ggplot2
pca <- plotPCA(vsd, intgroup = c("MHC_II","CDK12_status"), returnData = T)

#----------------------------------------------------------------------------------------
#------------Summarazing results and extracting significant gene lists-------------------
###---Results names:Factors compared in the design formula---
resultsNames(dds_DESeq) #Output: names to use in the shrink transformation

i <-  T
while (i) {
  result_shrink_Hemi_vs_Intact <- lfcShrink(dds_DESeq, coef = "",
                                            type = "apeglm", parallel = T)
  result_shrink_Homo_vs_Intact <- lfcShrink(dds_DESeq, coef = "",
                                            type = "apeglm", parallel = T)
  
  print(result_shrink_Hemi_vs_Intact)
  print(summary(result_shrink_Hemi_vs_Intact, alpha = 0.05))
  print(summary(result_shrink_Hemi_vs_Intact, alpha = 0.01))
  
  print(result_shrink_Homo_vs_Intact)
  print(summary(result_shrink_Homo_vs_Intact, alpha = 0.05))
  print(summary(result_shrink_Homo_vs_Intact, alpha = 0.01))
  i <- F
} 

###---MA plot: Mean of normalized counts vs log2 FC for all genes tested---
plotMA(result_shrink_Hemi_vs_Intact, ylim = c(-2,2), xlim = c(1,3),
       main = "MA plot - result_shrink:CDK12_status_CDK12_Het_vs_CDK12_WT")
plotMA(result_shrink_Homo_vs_Intact, ylim = c(-4,4), xlim = c(1,1e6), 
       main = "MA plot - result_shrink:CDK12_status_CDK12_Hom_vs_CDK12_WT")
#plotMA(result_shrink_Dupl_vs_Intact, ylim = c(-4,4), xlim = c(1,1e6), 
       #main = "MA plot - result_shrink:PTEN_SCNA_status_coding_Dupl_vs_Intact")

###---Create a tibble of results and extract sig. DE genes---
###---Cutoffs---
padj.cutoff_0.05 <- 0.05
padj.cutoff_0.01 <- 0.01

i <- T
while (i) {
  result_shrink_Hemi_vs_Intact_tb <- result_shrink_Hemi_vs_Intact %>%
    as.data.frame() %>%
    rownames_to_column(var = "gene") %>%
    as.data.frame() %>%
    arrange(padj)
  
  result_shrink_Homo_vs_Intact_tb <- result_shrink_Homo_vs_Intact %>%
    as.data.frame() %>%
    rownames_to_column(var = "gene") %>%
    arrange(padj)
  #result_shrink_Dupl_vs_Intact_tb <- result_shrink_Dupl_vs_Intact %>%
  #data.frame() %>%
  #rownames_to_column(var = "gene") %>%
  #as_tibble()
  i <- F
}

###---Subset the tibble to keep only significant genes(p < 0.05 and p < 0.01)---
i <- T
while (i) {
  res_0.05_shrink_Hemi_vs_Intac_tb <- result_shrink_Hemi_vs_Intact_tb %>%
    dplyr::filter(padj < padj.cutoff_0.05)  %>%
    arrange(padj)
  
  res_0.01_shrink_Hemi_vs_Intac_tb <- result_shrink_Hemi_vs_Intact_tb %>%
    dplyr::filter(padj < padj.cutoff_0.01)  %>%
    arrange(padj)
  
  res_0.05_shrink_Homo_vs_Intact_tb <- result_shrink_Homo_vs_Intact_tb %>%
    dplyr::filter(padj < padj.cutoff_0.05)  %>%
    arrange(padj)
  res_0.01_shrink_Homo_vs_Intact_tb <- result_shrink_Homo_vs_Intact_tb %>%
    dplyr::filter(padj < padj.cutoff_0.01)  %>%
    arrange(padj)
  i <- F
}

#----------------------------------------------------------------------------------------
#------------Graphich representation:Visualizing the results----------------------------- 
###---Generate the annotation table (tx2gene) using the ensembl ID to use in Graphs------
i <- T
while (i) {
  tx2gene <- annotables::grch38
  grch38annot <- tx2gene %>% 
    dplyr::select(ensgene, symbol, description, entrez, biotype) %>% 
    dplyr::distinct()
  i <- F
}

###---Generate the normalized gene column with matched ensembl ID---
##DESeq2 creates a matrix when you use the counts() function. First convert... 
##...normalized_counts to a data frame and transfer the row names to a new column... 
##...called "gene"
i <- T
while (i) {
  normalized_counts <- counts(dds_DESeq, normalized = T) %>%
    data.frame() %>%
    rownames_to_column(var = "gene")
  
    normalized_counts_tb <- normalized_counts %>% 
      left_join(grch38annot, by = c("gene" = "ensgene")) %>%
      distinct(gene, .keep_all = T)
    write.table(normalized_counts_tb, "", 
                row.names = T, 
                col.names = T, sep = "\t")
  i <- F
}

norm <- normalized_counts_tb
norm[,1] <- NULL

norm <- norm %>% 
  distinct(symbol, .keep_all = T) %>%
  drop_na(symbol)
rownames(norm) <- norm$symbol
norm <- norm %>% rownames_to_column(var = "gene")
norm[ , c(12:15)] <- NULL

write_tsv(norm, '', 
          col_names = T)

###---Generating the DE gene tables (p < 0.05 and p < 0.01)---
#This will bring in a column of gene symbols and write the results of DE sig. genes
i <- T
  while (i) {
  res_0.05_shrink_Hemi_vs_Intac_tb <- res_0.05_shrink_Hemi_vs_Intac_tb %>%
    left_join(grch38annot, by = c("gene"= "ensgene")) %>%
    distinct(gene, .keep_all = T)
  
  res_0.01_shrink_Hemi_vs_Intac_tb <- res_0.01_shrink_Hemi_vs_Intac_tb %>%
    left_join(grch38annot, by = c("gene" = "ensgene")) %>%
    distinct(gene, .keep_all = T) 
  
  #Write in the current directory
  write.table(res_0.05_shrink_Hemi_vs_Intac_tb, "", 
              row.names = T, col.names = T, sep = "\t")
  write.table(res_0.01_shrink_Hemi_vs_Intac_tb, "", 
              row.names = T, col.names = T, sep = "\t" )
  
  res_0.05_shrink_Homo_vs_Intact_tb <- res_0.05_shrink_Homo_vs_Intact_tb %>%
    left_join(grch38annot, by = c("gene" = "ensgene")) %>%
    distinct(gene, .keep_all = T) 
  
  res_0.01_shrink_Homo_vs_Intact_tb <- res_0.01_shrink_Homo_vs_Intact_tb %>%
    left_join(grch38annot, by = c("gene" = "ensgene")) %>%
    distinct(gene, .keep_all = T)
  
  #Write in the current directory
  write.table(res_0.05_shrink_Homo_vs_Intact_tb, "", 
              row.names = T, col.names = T, sep = "\t" )
  write.table(res_0.01_shrink_Homo_vs_Intact_tb, "", 
              row.names = T, col.names = T, sep = "\t" )
  i <- F
}

###---Extract normalized expression for significant genes---
i <- T
while (i) {
  vsd_mat_tb <- vsd_mat %>% as.data.frame() %>%
    rownames_to_column(var = "gene") %>%
    left_join(grch38annot, by = c("gene" = "ensgene")) %>%
    distinct(gene, .keep_all = T)
  
  vsd_mat_Hem <- vsd_mat_tb %>% 
    dplyr::filter(gene %in% res_0.05_shrink_Hemi_vs_Intac_tb$gene) %>%
    distinct(symbol, .keep_all = T) %>%
    left_join(res_0.05_shrink_Hemi_vs_Intac_tb, by = c("gene" = "gene")) %>%
    dplyr::arrange(padj) %>%
    drop_na(symbol.x)
  
  #name row with gene symbol(HUGO)
  rownames(vsd_mat_Hem) <- vsd_mat_Hem$symbol.x
  
  vsd_mat_Hom <- vsd_mat_tb %>% 
    dplyr::filter(gene %in% res_0.05_shrink_Homo_vs_Intact_tb$gene) %>%
    distinct(symbol, .keep_all = T) %>%
    left_join(res_0.05_shrink_Homo_vs_Intact_tb, by = c("gene" = "gene")) %>%
    dplyr::arrange(padj)
  
  #name row with gene symbol(HUGO)
  rownames(vsd_mat_Hom) <- vsd_mat_Hom$symbol.x
  
  vsd_mat_total <- vsd_mat_tb %>%
    distinct(symbol, .keep_all = T) %>%
    drop_na() 
  rownames(vsd_mat_total) <- vsd_mat_total$symbol
  i <- F
}

#---Plotting sig. DE genes:Using 'plotCounts()' function---------------------------------
##"plotCounts()" function requires exactly original input to DESeq2(Ensembl IDs)

HLA_A$ <- factor(HLA_A$, levels = c("High", "Low"))
###---plot of counts(DE transcript)---
ggplot(HLA_A, aes(x = , y = count))+ 
  ggtitle("plotCounts: STAT4")    +
  scale_fill_manual(values = c("red2", "dodgerblue2"))  +
  geom_point(position = position_jitter(w = 0.4,h = 0), alpha=0.3)+
  theme_classic()                                       +
  geom_boxplot(aes(fill = ), outlier.size = 0.07, alpha = 0.7, size = 0.5)                            +
  theme(axis.text = element_text(size = 10, face = 'bold'),
        axis.title = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        plot.title = element_text(hjust = 0, size = 10),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13))          +
  scale_y_log10() + stat_compare_means(aes(group = ), paired = F, method = "t.test", 
                                       label = "p.signif", hide.ns = T, ref.group = "Low") +
  geom_text_repel(aes(label = a), max.overlaps = 5)

a <- row.names(HLA_A)

###---Heatmap of the first 50 DE genes---
##Set up the color pallet
heat_colors <- maPalette(low = "blue", high = "red", mid = 'white', k = 1000) #Valor de K == 50 

##Run the "pheatmap" using the metadata data frame for the annotation
##Use "count()" to make a index of metadata rows to map the samples/groups in the heatmap
meta <- dplyr::select(metadata, MHC_II) %>%
  mutate(count = row_number(MHC_II))
meta2 <- dplyr::select(metadata, c(MHC_II, CDK12_status, Prior.Treatment, Chemo.Regimen.Category,
                                   Tumor.Site))

##Run the "pheatmap()"
pheatmap(vsd_mat_Hem[1:50, 2:11], 
         color         = heat_colors,
         cluster_rows  = T,
         cluster_cols  = F, 
         show_rownames = T,
         show_colnames = T,
         annotation    = meta2, 
         border_color  = NA, 
         fontsize      = 10, 
         scale         = "row", 
         fontsize_row  = 12, 
         height        = 20,
         main          = "Heatmap: mCRPC_CDK12_MHC_I High vs Low")#,
         annotation_colors = list(CDK12 = c(CDK12_HET_LOSS = "Green", CDK12_WT = "lightpink")))


#if cluster_cols == T, save as 900x900, else(cluster == F), save as 800x800
pheatmap(norm_HOMO_vc_Intact_sig[1:100, c(43:489)], 
         color         = heat_colors,
         cluster_rows  = T,
         cluster_cols  = F, 
         show_rownames = T,
         show_colnames = F,
         annotation    = meta2, 
         border_color  = NA, 
         fontsize      = 10, 
         scale         = "row", 
         fontsize_row  = 12, 
         height        = 20,
         main          = "Heatmap: TCGA_PRAD-CDK12_Hom vs WT",
         annotation_colors = list(CDK12 = c(CDK12_HOM_LOSS = "Blue", CDK12_WT = "lightpink")))

###---Heatmap of a subset of the DE genes---
##The subset of genes can be set up by one's interest using the gene symbols in the... 
##...'norm_counts_tb_sig'
##Exemple
hla <- c("HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DQB2", 
         "HLA-DRA", "HLA-DRB1", "HLA-DRB5", "HLA-DRB6", "HLA-DPA1", "HLA-DPB1", "HLA-DMA", 
         "HLA-DMB")
c("HLA-A", "HLA-B", "HLA-C")#
chem_receptor <- c("CCL1", "CCL2", "CCL7", "CCL8", "CCL11", "CCL13", "CCL17", 
                   "CCL19", "CCL20", "CCL21", "CCL22", "CCL24", "CCL25", "CCL26", 
                   "CXCL1", "CXCL2", "CXCL3", "CXCL5", "CXCL6", "CXCL8", "CXCL9", 
                   "CXCL10", "CXCL11", "CXCL12", "CXCL13", "CXCL14", "CXCL16")
chem_receptor <- c("CCL1", "CCL2", "CCL7", "CCL8", "CCL11", "CCL13", "CCL17", 
                   "CCL19", "CCL20", "CCL21", "CCL22", "CCL24", "CCL25", "CCL26", 
                   "CXCL1", "CXCL2", "CXCL3", "CXCL5", "CXCL6", "CXCL8", "CXCL9", 
                   "CXCL10", "CXCL11", "CXCL12", "CXCL14", "CXCL16")

icp <- c("LAG3", "IDO1", "PDCD1", "CD274", "IFNG", "HAVCR2")

pheatmap(vsd_mat_total[hla, 2:49], 
         color         = heat_colors,
         cluster_rows  = T,
         cluster_cols  = F, 
         show_rownames = T,
         show_colnames = F,
         annotation    = meta2, 
         border_color  = NA, 
         fontsize      = 15, 
         scale         = "row", 
         fontsize_row  = 15, 
         height        = 20,
         main          = "Heatmap: PRAD_MHC_II")#,
         annotation_colors = list(CDK12 = c(CDK12_HET_LOSS = "Green", CDK12_HOM_LOSS = "Blue")))


#----------------------------------------------------------------------------------------
#------------Functional Analysis:clusterProfiler----------------------------------------- 
###---Over-Representation(OR) method:Hypergeometric test---
i <- T
while (i) {
  ##Create background dataset for hypergeometric testing using all genes tested for significance in the results                 
  all_genes <- as.character(result_shrink_Hemi_vs_Intact_tb$gene)
  #all_genes2 <- as.character(result_shrink_Homo_vs_Intact_tb$gene)
  
  ##Extract significant results
  sig_up <- res_0.05_shrink_Hemi_vs_Intac_tb %>%
    dplyr::filter(res_0.05_shrink_Hemi_vs_Intac_tb$log2FoldChange > 0)
  sig_down <- res_0.05_shrink_Hemi_vs_Intac_tb %>%
    dplyr::filter(res_0.05_shrink_Hemi_vs_Intac_tb$log2FoldChange < 0)
  sig_genes_up <- as.character(sig_up$gene)
  sig_genes_down <- as.character(sig_down$gene)
  
  ##Running OR enrichment analyses 
  ego_up <- enrichGO(gene  = sig_genes_up, 
                     universe = all_genes,
                     keyType  = "ENSEMBL",
                     OrgDb    = org.Hs.eg.db, 
                     ont      = "BP", #Go terms:Biological Process(BP), Cellular Component(CC), Molecular Function(MF), or (ALL)
                     pAdjustMethod = "BH", 
                     qvalueCutoff  = 0.05, 
                     readable = TRUE) #Gives the gene symbol as output(matched with genes IDs)
  
  ego_down <- enrichGO(gene = sig_genes_down, 
                       universe = all_genes,
                       keyType  = "ENSEMBL",
                       OrgDb    = org.Hs.eg.db, 
                       ont      = "BP", 
                       pAdjustMethod = "BH", 
                       qvalueCutoff  = 0.05, 
                       readable      = TRUE)
  
  ##use 'simplify()' to remove redundant terms
  ego_up_simply <- simplify(ego_up, cutoff   = 0.7, by = "p.adjust", select_fun = min)
  ego_down_simply <- simplify(ego_down, cutoff = 0.7, by = "p.adjust", select_fun = min)
  #ego_simply_new_up <- mutate(ego_up_simply, 
                           #richFactor  = Count / as.numeric(sub("/\\d+", "", BgRatio)))
  #ego2_simply_new_down <- mutate(ego_down_simply, 
                            #richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
  i <- F
}

##Check out the enriched results descriptions
ego_up_simply@result$Description
ego_down_simply@result$Description

###---Visualizing the OR clusterProfiler results---
##dotplot()
dotplot(ego_up_simply, showCategory  = 40, x  = 'richFactor', font.size = 15, title = "enrich_TCGA_PRAD_CDK12_Het_WT_UP")
dotplot(ego2_simply_new_down, showCategory = 40, x  = 'richFactor', font.size = 15, title = "enrich_TCGA_PRAD_CDK12_Het_DOWN")

##enrichmentp GO plot:Relationship between the top most sig. enrichement GO terms (padj.)
##Add similarity matrix to the termsim slot of enrichment result
ego_simply_up <- enrichplot::pairwise_termsim(ego_up_simply)
ego2_simply_down <- enrichplot::pairwise_termsim(ego_down_simply)
##Enrichmap cluster
emapplot(ego_simply_up, showCategory  = 50, color  = 'p.adjust')
emapplot(ego2_simply_down, showCategory = 50,  color = 'p.adjust')

##category plot:Relationship of the top five most sig. GO terms and fold changes associeted with it 
##To color genes by log2 fold changes, we need to extract the log2 fold changes from our results table creating a named vector
sig_foldchanges_up <- sig_up$log2FoldChange
sig_foldchanges_down <- sig_down$log2FoldChange
names(sig_foldchanges_up) <- sig_up$gene
names(sig_foldchanges_down) <- sig_down$gene
##If some of the high fold changes are getting drowned out due to a large range, you could set a maximum fold change value
sig_foldchanges_up <- ifelse(sig_foldchanges_up >  2, 2, sig_foldchanges_up)
sig_foldchanges_up <- ifelse(sig_foldchanges_up < -2, -2, sig_foldchanges_up)

sig_foldchanges_down <- ifelse(sig_foldchanges_down > 2, 2, sig_foldchanges_down)
sig_foldchanges_down <- ifelse(sig_foldchanges_down < -2, -2, sig_foldchanges_down)

##cnetplot() - by default gives the top 5 sig. terms (padj.)
cnetplot(ego_simply_up, 
         categorySize = "pvalue", 
         showCategory = 5, 
         foldChange   = sig_foldchanges_up, 
         vertex.label.font = 6, node_label = 'all')
cnetplot(ego2_simply_down, 
         categorySize = "pvalue", 
         showCategory = 2, 
         foldChange   = sig_foldchanges_down, 
         vertex.label.font = 6, node_label = 'all')

#Description and subset
ego_simply_sub <- ego_simply_up
ego_simply_sub@result$Description
ego_simply_sub@result <- ego_simply_sub@result[c(3), ]#personal criteria
cnetplot(ego_simply_sub, 
         categorySize = "pvalue", 
         showCategory = 1, 
         foldChange   = sig_foldchanges_up, 
         vertex.label.font = 6, node_label = 'all')

###---Functional class scoring / GSEA method---
##gse_GO():GSEA method
#Merge the annotation table (tx2gene) dataframe with the results(all transcripts)
i <- T
while (i) {
  res_ids <- left_join(result_shrink_Hemi_vs_Intact_tb, grch38annot, by  = c("gene" = "ensgene"))
  res2_ids <- left_join(result_shrink_Homo_vs_Intact_tb, grch38annot, by = c("gene" = "ensgene"))
  #Remove any NA values
  res_entrez <- dplyr::filter(res_ids, entrez   != "NA")
  res2_entrez <- dplyr::filter(res2_ids, entrez != "NA")
  #Remove entrez duplicated
  res_entrez <- res_entrez[which(duplicated(res_entrez$entrez)  == F), ]
  res2_entrez <- res2_entrez[which(duplicated(res2_entrez$entrez) == F), ]
  #Extract foldchanges
  foldchanges <- res_entrez$log2FoldChange
  foldchanges2 <- res2_entrez$log2FoldChange
  #Name each fold change with the corresponding Entrez ID
  names(foldchanges) <- res_entrez$entrez
  names(foldchanges2) <- res2_entrez$entrez
  #sort the foldchanges in decreasing order
  foldchanges <- sort(foldchanges, decreasing   = TRUE)
  foldchanges2 <- sort(foldchanges2, decreasing = TRUE)
  
  ##GSEA using GO terms
  gsego <- gseGO(geneList     = foldchanges, 
                 OrgDb        = org.Hs.eg.db, 
                 ont          = 'MF', 
                 nPermSimple  = 100000, 
                 minGSSize    = 10, 
                 pvalueCutoff = 0.05,
                 verbose = FALSE, 
                 seed    = T, 
                 eps     = 1e-10) 
  gsego2 <- gseGO(geneList    = foldchanges2, 
                  OrgDb        = org.Hs.eg.db, 
                  ont          = 'MF', 
                  nPermSimple  = 100000, 
                  minGSSize    = 10, 
                  pvalueCutoff = 0.05,
                  verbose = FALSE) 
  ##use 'simplify()' to remove redundant terms
  gsego_simply <- simplify(gsego, cutoff   = 0.7, by = "p.adjust", select_fun = min)
  gsego_simply2 <- simplify(gsego2, cutoff = 0.7, by = "p.adjust", select_fun = min)
  i <- F
}

##Check out the enriched results descriptions
gsego_simply@result$Description
gsego_simply2@result$Description

##Write in the current directory
write.table(gsego_simply@result, "Res/GSEA_GO_MF_Condition_IC_vs_Control.txt", 
            row.names = TRUE, col.names = TRUE, sep = "\t" )
write.table(gsego_simply2@result, "Res/GSEA_GO__MF_Condition_Dx.IONP_vs_Control.txt", 
            row.names = TRUE, col.names = TRUE, sep = "\t" )

##plot the gse_GO() enriched results
##slice the list between two groups(NES + / -)
ewp <- arrange(gsego_simply, desc(abs(NES))) %>%
  group_by(sign(NES)) %>%
  slice(1:35)

ggplot(ewp, showCategory = 120,
       aes(NES, fct_reorder(Description, NES), fill = p.adjust)) +
  geom_col() +
  scale_fill_gradientn(colours = c('#b3eebe', "#46bac2", '#371ea3'),
                       guide = guide_colorbar(reverse = T)) +
  #theme_dose(12) +
  theme_bw()+
  theme(plot.title   = element_text( size = 20, hjust = 0.5),
        axis.text    = element_text(size = 20, color = 'black'),
        axis.title.x = element_text(size = 15),
        legend.title = element_text(size = 20),
        legend.text  = element_text(size = 15))+
  xlab("Normalized Enrichment Score") +
  ylab(NULL) +
  ggtitle("GSEA_GO(MF)_Condition_IC_vs_Control")

ewp2 <- arrange(gsego_simply2, desc(abs(NES))) %>%
  group_by(sign(NES)) %>%
  slice(1:26)

ggplot(ewp2, showCategory=120,
       aes(NES, fct_reorder(Description, NES), fill = p.adjust)) +
  geom_col() +
  scale_fill_gradientn(colours = c('#b3eebe', "#46bac2", '#371ea3'),
                       guide = guide_colorbar(reverse = T)) +
  theme_bw()+
  theme(plot.title   = element_text(size = 20, hjust = 0.5),
        axis.text    = element_text(size = 20, color = 'black'),
        axis.title.x = element_text(size = 15),
        legend.title = element_text(size = 20),
        legend.text  = element_text(size = 15))+
  xlab("Normalized Enrichment Score") +
  ylab(NULL) +
  ggtitle("GSEA_GO(MF)_Condition_Dx.IONP_vs_Control")

##GSEA():GSEA method with Molecular Signature Data Base (MSigDB) gene set
##MSigDB collections
msigdbr <- msigdbr_collections()

##Set up the deseired gene set collections 
i <- T
while (i) {
  m_t2g_H <- msigdbr(species = 'Homo sapiens', category = 'H') %>%
    dplyr::select(gs_name, entrez_gene)
  m_t2g_C2_REACTOME <- msigdbr(species = 'Homo sapiens', category = 'C2', subcategory = 'CP:REACTOME') %>%
    dplyr::select(gs_name, entrez_gene)
  m_t2g_C5 <- msigdbr(species = 'Homo sapiens', category = 'C5', subcategory = 'GO:BP') %>%
    dplyr::select(gs_name, entrez_gene)
  m_t2g_C6 <- msigdbr(species = 'Homo sapiens', category = 'C6') %>%
    dplyr::select(gs_name, entrez_gene)
  m_t2g_C7 <- msigdbr(species = 'Homo sapiens', category = 'C7', subcategory = 'IMMUNESIGDB') %>%
    dplyr::select(gs_name, entrez_gene)
  m_t2g_C8 <- msigdbr(species = 'Homo sapiens', category = 'C8') %>%
    dplyr::select(gs_name, entrez_gene)
  i <- F
}
##GSEA()
i <- T
while (i) {
  gsea_Hemi <- GSEA(foldchanges, exponent = 1,
                    TERM2GENE    = m_t2g_C2_REACTOME,
                    nPermSimple  = 100000,
                    verbose      = T,
                    pvalueCutoff = 0.05,
                    by = 'fgsea', eps = 1e-10)
  gsea_Hom <- GSEA(foldchanges2, exponent = 1,
                   TERM2GENE    = m_t2g_C7,
                   nPermSimple  = 100000,
                   verbose      = T,
                   pvalueCutoff = 0.05,
                   by = 'fgsea', eps = 1e-10)
  
  ##Transform the GSEA() results as data.frame
  res_gsea_Hemi <- as.data.frame(gsea_Hemi)
  res_gsea_Homo <- as.data.frame(gsea_Hom)
  i <- F
}

##plot the GSEA() enriched results
##slice the list between two groups(NES + / -)
ewp <- arrange(res_gsea_Hemi, desc(abs(NES))) %>%
  dplyr::group_by(sign(NES)) %>%
  slice(1:30)
ggplot(ewp, showCategory = 28, aes(NES, fct_reorder(Description, NES), 
                                   fill=p.adjust))                     +
  geom_col()      +
  scale_fill_gradientn(colours = c('#b3eebe', "#46bac2", '#371ea3'),
                       guide = guide_colorbar(reverse = T))            +
  #theme_dose(12) +
  theme_bw()      +
  theme(plot.title   = element_text(size = 20, hjust = 0.5),
        axis.text    = element_text(size = 15, color = 'black'),
        axis.title.x = element_text(size = 15),
        legend.title = element_text(size = 20),
        legend.text  = element_text(size = 15)) +
  xlab("Normalized Enrichment Score")           +
  ylab(NULL)                                    +
  ggtitle("Reactome: CDK12_MHC_II High vs Low")

ewp2 <- arrange(res_gsea_Homo, desc(abs(NES))) %>%
  dplyr::group_by(sign(NES)) %>%
  slice(1:30)
ggplot(ewp2, showCategory = 60, aes(NES, fct_reorder(Description, NES), 
                                    fill = p.adjust))                  +
  geom_col()      +
  scale_fill_gradientn(colours = c('#b3eebe', "#46bac2", '#371ea3'),
                       guide = guide_colorbar(reverse = T))            +
  #theme_dose(12) +
  theme_bw()      +
  theme(plot.title   = element_text(size = 20, hjust = 0.5),
        axis.text    = element_text(size = 15, color = 'black'),
        axis.title.x = element_text(size = 15),
        legend.title = element_text(size = 20),
        legend.text  = element_text(size = 15)) +
  xlab("Normalized Enrichment Score")           +
  ylab(NULL)                                    +
  ggtitle("Reactome Condition_Dx.IONP_vs_Control")