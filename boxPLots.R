#adrenal 
acc <- read_tsv("TCGA_adrenal_gland/results/expression_Adrenal.txt",
                col_names = T)
acc <- acc %>%
  dplyr::mutate(Tumor_type = str_replace_all(acc$Tumor_type, pattern = "TCGA_Adrenal", replacement = "ACC"))

#Bladder
blca <- read_tsv("TCGA_bladder/results/expression_Bladder.txt",
                col_names = T)
blca <- blca %>%
  dplyr::mutate(Tumor_type = str_replace_all(blca$Tumor_type, pattern = "TCGA_Adrenal", replacement = "BLCA"))

##Brain
#GBM
gbm <- read_tsv("TCGA_brain/Results/results_GBM/expression_GBM.txt",
                col_names = T)
gbm <- read_tsv("TCGA_brain/Results/brain/expression_Brain.txt",
                col_names = T)
gbm <- gbm %>%
  dplyr::mutate(Tumor_type = str_replace_all(gbm$Tumor_type, pattern = "TCGA_GBM", replacement = "GBM"))

lgg <- read_tsv("TCGA_brain/Results/results_LGG/expression_LGG.txt",
                col_names = T)

#breast 
brca <- read_tsv("TCGA_breast/results/expression_Breast.txt",
                 col_names = T)

#cervix
cervix <- read_tsv("TCGA_cervix/results/expression_cervix.txt",
                   col_names = T)
#colorectal
colo <- read_tsv("TCGA_colorectal/results_COAD/expression_COAD.txt",
                 col_names = T)
colo2 <- read_tsv("TCGA_colorectal/results_READ/expression_READ.txt",
                  col_names = T)
#esophagus
eso <- read_tsv("TCGA_esophagus/results/expression_esophagus.txt",
                col_names = T)
# eye
eye <- read_tsv("TCGA_eye/results/expression_eye.txt",
                col_names = T)

#heead and neck 
hand <- read_tsv("TCGA_head_and_neck/results/expression_head_and_neck.txt",
                 col_names = T)

#kidney
kidney <- read_tsv("TCGA_kidney/results_KICH/expression_KICH.txt",
                   col_names = T)
kidney2 <- read_tsv("TCGA_kidney/results_KIRC/expression_KIRC.txt",
                    col_names = T)
kidney3 <- read_tsv("TCGA_kidney/results_KIRP/expression_KIRP.txt",
                    col_names = T)
#liver
liver <- read_tsv("TCGA_liver/results/expression_liver.txt",
                    col_names = T)

#lung
LUAD <- read_tsv("TCGA_lung/results_LUAD/expression_LUAD.txt",
                    col_names = T)
lusc <- read_tsv("TCGA_lung/results_LUSC/expression_LUSC.txt",
                    col_names = T)

#lymoh node 
lymph <- read_tsv("TCGA_lymph_nodes/results/expression.txt",
                 col_names = T)

#ovary 
ov <- read_tsv("TCGA_ovary/results/expression.txt",
                  col_names = T)
#pancreas
pancreas <- read_tsv("TCGA_pancreas/results/expression.txt",
                  col_names = T)

#pleura 
pleura <- read_tsv("TCGA_pleura/results/expression.txt",
                  col_names = T)
#prostate 
prad <- read_tsv("TCGA_prostate/results/expression.txt",
                  col_names = T)
#skin 
skin <- read_tsv("TCGA_skin/results/expression.txt",
                  col_names = T)
#soft  
soft <- read_tsv("TCGA_soft_tissue/results/expression.txt",
                  col_names = T)
#stomach 
stom <- read_tsv("TCGA_stomach/results/expression_STAD.txt",
                  col_names = T)

#testis 
testis <- read_tsv("TCGA_testis/results/expression_TGCT.txt",
                  col_names = T)

#thymus 
thymus <- read_tsv("TCGA_thymus/results/expression_THYM.txt",
                  col_names = T)
#thyroid 
thy <- read_tsv("TCGA_thyroid/results/expression_THCA.txt",
                  col_names = T)
#utherus
uth1 <- read_tsv("TCGA_uterus/results_UCEC/expression_UCEC.txt",
                  col_names = T)
uth2 <- read_tsv("TCGA_uterus/results_UCS/expression_UCS.txt",
                 col_names = T)

table <- bind_rows(acc, blca, brca, cervix, colo, colo2, eso, eye, gbm,
                   hand, kidney, kidney2, kidney3, lgg, liver, LUAD, lusc,
                   lymph, ov, pancreas, pleura, prad, skin, soft, stom,
                   testis, thy, thymus, uth1, uth2)
table <- bind_rows(acc, blca, brca, cervix, colo, colo2, eso, eye, gbm,
                   hand, kidney, kidney2, kidney3, liver, LUAD, lusc,
                   lymph, ov, pancreas, pleura, prad, skin, soft, stom,
                   testis, thy, thymus, uth1, uth2)

names <- read_tsv("../../../Downloads/mmc2.txt", col_names = T) 
names <- names[,c(1,2)]
names <- names %>%
  dplyr::mutate(ID = str_replace_all(names$ID, pattern = "[-]", replacement = "."))

table <- table %>%
  left_join(names, by = c("ID" = "ID"))

table[,4] <- NULL
table <- table %>%
  dplyr::mutate(PTEN_SCNA_status = str_replace(table$PTEN_SCNA_status, "Hemi", "HemDel"))

table <- table %>%
  dplyr::mutate(PTEN_SCNA_status = str_replace(table$PTEN_SCNA_status, "Homo", "HomDel"))

table_teste <- table %>%
  dplyr::mutate(TCGA_Study = str_replace(table$TCGA_Study, "LGG", "Brain"))
table_teste <- table %>%
  dplyr::mutate(TCGA_Study = str_replace(table$TCGA_Study, "GBM", "Brain"))

table3 <- table %>%
  drop_na() %>%
  dplyr::filter(PTEN_SCNA_status == c("HemDel", "HomDel", "Intact")) %>%
  dplyr::filter( TCGA_Study != "UVM", TCGA_Study != "DLBC", TCGA_Study != "KICH",
                 TCGA_Study != "KIRP", TCGA_Study != "PAAD", TCGA_Study != "PCPG",
                 TCGA_Study != "READ", TCGA_Study != "TGCT", TCGA_Study != "THCA", 
                 TCGA_Study != "THYM", TCGA_Study != "UCEC", TCGA_Study != "ACC",
                 TCGA_Study != "BLCA", TCGA_Study != "CESC", TCGA_Study != "COAD",
                 TCGA_Study != "ESCA", TCGA_Study != "KIRC", TCGA_Study != "LGG",
                 TCGA_Study != "LUAD", TCGA_Study != "LUSC", TCGA_Study != "MESO",
                 TCGA_Study != "SKCM") %>%
  dplyr::filter(Gene != "IFNG", Gene != "PDCD1", Gene != "PTEN", Gene != "HAVCR2")

table3$PTEN_SCNA_status <- factor(table3$PTEN_SCNA_status, levels = c("HomDel","HemDel","Intact"))

table4 <- table %>%
  drop_na() %>%
  dplyr::filter(PTEN_SCNA_status == c("HemDel", "HomDel", "Intact")) %>%
  dplyr::filter(TCGA_Study != "UVM", TCGA_Study != "DLBC", TCGA_Study != "KICH",
                TCGA_Study != "KIRP", TCGA_Study != "PAAD", TCGA_Study != "PCPG",
                TCGA_Study != "READ", TCGA_Study != "TGCT", TCGA_Study != "THCA", 
                TCGA_Study != "THYM", TCGA_Study != "UCEC") %>%
  dplyr::filter(Gene == "PTEN") 

table6$PTEN_SCNA_status <- factor(table6$PTEN_SCNA_status, levels = c("HomDel", "HemDel", "Intact"))

color <- c("dodgerblue2","red2","Darkgoldenrod2","forestgreen","lightblue")
library(egg)
library(ggplot2)
library(ggsignif)
library(ggpubr)


table6$TCGA_Study <- factor(table6$TCGA_Study, 
                            levels = c("SKCM", "LUSC", "SARC", 
                                       "TGCT", "BLCA", "OV",
                                       "Brain", "ESCA", "PRAD",
                                       "READ", "BRCA", "CESC",
                                       "LUAD", "STAD", "HNSC",
                                       "LIHC", "COAD", "KIRC",
                                       "PAAD", "KIRP", "THYM",
                                       "THCA", "PCPG"))

#Look last code


library(ggpubr)
library(rstatix)
stat.test
symnum.args <- list(cutpoints = c(0, 0.05, 1), symbols = c( "*", "ns"))

legend.position = c(0.9, 0.2),
legend.key.height = unit(2, "mm")
legend.position = "none"
p + 
  geom_boxplot(aes(fill = PTEN_SCNA_status), outlier.size = 0.07, size = 0.5) + 
  facet_grid(rows = vars(table3$Gene),cols = vars(table3$TCGA_Study), switch = "x",
            scales = "free") +
  ggtitle("") +
  xlab("") +
  ylab("Transcript counts (Log2)") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 0),
        axis.text.y = element_text(size = 5),
        legend.position = "none",
        strip.text = element_text(size = 6),
        strip.text.x = element_text(angle = 0),
        strip.placement = "bottom",
        strip.background.x=element_rect(color = "white",  fill=NULL), 
        strip.background.y=element_rect(color = NA)) +
  scale_fill_manual(values = color) +
  stat_compare_means(aes(group = PTEN_SCNA_status), 
                     paired = F, method = "t.test", 
                     label = "p.signif", hide.ns = T, ref.group = "Intact")



table3 <- table_teste %>%
  dplyr::filter(Gene == c("PTEN"))

table3 <- table3 %>%
  drop_na()

table_teste_2 <- table_teste %>%
  dplyr::filter(Gene == c("PTEN"))

table_teste_2 <- table_teste_2 %>%
  dplyr::mutate(TCGA_Study = str_replace(table_teste_2$TCGA_Study, "LGG", "Brain"))
table_teste_2 <- table_teste_2 %>%
  dplyr::mutate(TCGA_Study = str_replace(table_teste_2$TCGA_Study, "GBM", "Brain"))







library(tidyverse)
library(egg)
table5 <- read_tsv("Downloads/table_figure1.txt", col_names = T)

table6 <- table5 %>%
  drop_na(TCGA_Study) %>%
  filter(TCGA_Study != "ACC", TCGA_Study != "MESO", TCGA_Study != "DLBC",
         TCGA_Study != "UCEC", TCGA_Study != "UCS", TCGA_Study != "UVM",
         TCGA_Study != "NA") %>%
  filter(Gene == "PTEN")

table6$PTEN_SCNA_status <- factor(table6$PTEN_SCNA_status, 
                                  levels = c("HemDel","HomDel","Dupl","Intact"))
table6$TCGA_Study <- factor(table6$TCGA_Study, 
                            levels = c("PCPG","THCA", "THYM",
                                       "KIRP", "PAAD", "KIRC",
                                       "COAD", "LIHC", "HNSC",
                                       "STAD", "LUAD", "CESC",
                                       "BRCA", "READ", "PRAD",
                                       "ESCA", "Brain", "OV",
                                       "BLCA", "TGCT", "SARC",
                                       "LUSC", "SKCM"))


table6 <- table6 %>%
  drop_na(TCGA_Study)

ggplot(table6)+
  geom_bar(mapping = aes(y=TCGA_Study, fill = PTEN_SCNA_status), position = "fill")+
  scale_fill_manual(values=c("dodgerblue2","red2","Darkgoldenrod2","forestgreen","lightblue"))+
  #facet_wrap(~ EC , nrow = 2)+
  theme_article()+
  #scale_x_discrete('')+
  scale_x_continuous(name='Percentage' ,labels = scales::percent)+
  theme(axis.title=element_text(size=10),
        axis.text.x = element_text(),
        legend.position = "bottom",
        legend.title = element_text(size = 8)) +
  xlab("") + 
  ylab("")

write_tsv(table6, "Desktop/PTEN/PAPER/table_figure1.txt", col_names = T)











#Figure 3a
p <- ggplot(table6, mapping = aes(x = PTEN_SCNA_status, y = Count)) 
p + 
  geom_boxplot(aes(fill = PTEN_SCNA_status), outlier.size = 0.07, size = 0.5) + 
  #facet_grid(cols = vars(table6$TCGA_Study), switch = "x") +
  ggtitle("") +
  xlab("") +
  ylab("Transcript counts (Log2)") +
  theme_article() +
  theme(axis.text.x = element_text(size = 0),
        axis.text.y = element_text(size = 7, color = "Black"),
        #legend.position = c(0.9, 0.1),
        legend.position = "bottom",
        #legend.key.height = unit(2, "mm"),
        legend.title = element_text(size = 9),
        legend.box.background = element_rect(color = "white"),
        strip.text = element_text(size = 6),
        strip.text.x = element_text(angle = 0),
        strip.placement = "bottom", 
        strip.background.x=element_rect(color = NA,  fill=NA), 
        strip.background.y=element_rect(color = "black",  fill=NA)) +
  scale_fill_manual(values = color) +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", y.position = c(22, 23, 24),
                     tip.length = 0.01, bracket.shorten = .04)
  

stat.test <- table6 %>%
  t_test(Count ~ PTEN_SCNA_status, ref.group = "Intact") %>%
  add_significance()
stat.test <- stat.test %>% add_xy_position(x = "PTEN_SCNA_status") 

stat.test$y.position


table6_HemDel <- table6 %>% filter(PTEN_SCNA_status == "HemDel")
table6_HomDel <- table6 %>%  filter(PTEN_SCNA_status == "HomDel")
table6_Depl <- table6 %>%  filter(PTEN_SCNA_status == "Dupl")
table6_Intact <- table6 %>%  filter(PTEN_SCNA_status == "Intact")
summary(table6_HemDel)
summary(table6_HomDel)
summary(table6_Depl)
summary(table6_Intact)





#figure 4

table <- read_tsv("Downloads/mmc2.txt")
table2 <- read_tsv("Desktop/PTEN/metadata/PTEN_groups_by_tumor.txt",
                   col_names = T)
table2 <- table2 %>%
  dplyr::mutate(ID = str_replace_all(table2$ID, pattern = "[.]", replacement = "-"))


table2 <- table2[,c("ID", "PTEN_SCNA_status")]
table <- table %>%
  left_join(table2)

table <- table %>%
  dplyr::filter(PTEN_SCNA_status != c("Dupl"))

table <- table %>%
  dplyr::filter(PTEN_SCNA_status != c("HomDel","HemDel", "NA"))

table$PTEN_SCNA_status <- factor(table$PTEN_SCNA_status, levels = c("HomDel","HemDel","Intact"))
table$PTEN_SCNA_status <- factor(table$PTEN_SCNA_status, levels = c("Dupl","Intact"))

library(ggplot2)
#plot
color <- c("red2", "dodgerblue2", "forestgreen") #Darkgoldenrod2
color <- c("Darkgoldenrod2", "forestgreen") 

b4 <- ggplot(table, mapping = aes(x = PTEN_SCNA_status, 
                                  y = Intratumor_Heterogeneity, 
                                  fill = PTEN_SCNA_status )) 
ith <- b4 + geom_point(aes(color = PTEN_SCNA_status),
                       position = position_jitter(w = 0.3,h = 0), alpha = 0.3, ) +
  geom_boxplot(aes(fill = PTEN_SCNA_status),
               outlier.size = 0.07, size = 0.5) +
  xlab("") +
  ylab("ITH levels") +
  theme_article() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 7, color = "Black"),
        #legend.position = c(0.9, 0.1),
        legend.position = "none",
        #legend.key.height = unit(2, "mm"),
        legend.title = element_text(size = 0),
        legend.box.background = element_rect(color = "white"),
        strip.text = element_text(size = 6),
        strip.text.x = element_text(angle = 0),
        strip.placement = "bottom", 
        strip.background.x=element_rect(color = NA,  fill=NA), 
        strip.background.y=element_rect(color = "black",  fill=NA)) +
  scale_fill_manual(values = color) +
  scale_color_manual(values = color) +
  stat_compare_means(comparisons = list(c("HomDel", "HemDel"), 
                                        c("HemDel","Intact"), 
                                        c("HomDel","Intact")), 
                     label.y = c(0.25, 0.30, 0.15) + 
                       stat_compare_means(label.y = 2))

ith



a4 <- ggplot(table, mapping = aes(x = PTEN_SCNA_status, 
                                  y = log(Nonsilent_Mutation_Rate), 
                                  fill = PTEN_SCNA_status )) 

nsilMut <- a4 + geom_point(aes(color = PTEN_SCNA_status),
                           position = position_jitter(w = 0.3,h = 0), alpha = 0.3) +
  geom_boxplot(aes(fill = PTEN_SCNA_status),
               outlier.size = 0.07, size = 0.5) +
  xlab("") +
  ylab("Nonsilent Mutation") +
  theme_article() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 7, color = "Black"),
        #legend.position = c(0.9, 0.1),
        legend.position = "none",
        #legend.key.height = unit(2, "mm"),
        legend.title = element_text(size = 0),
        legend.box.background = element_rect(color = "white"),
        strip.text = element_text(size = 6),
        strip.text.x = element_text(angle = 0),
        strip.placement = "bottom", 
        strip.background.x=element_rect(color = NA,  fill=NA), 
        strip.background.y=element_rect(color = "black",  fill=NA)) +
  scale_fill_manual(values = color) +
  scale_color_manual(values = color) +
  stat_compare_means(comparisons = list(c("HomDel", "HemDel"), 
                                        c("HemDel","Intact"), 
                                        c("HomDel","Intact")), 
                     label.y = c(0.25, 0.30, 0.15) + 
                       stat_compare_means(label.y = 2))
nsilMut 



c4 <- ggplot(table, mapping = aes(x = PTEN_SCNA_status, 
                                  y = Homologous_Recombination_Defects, 
                                  fill = PTEN_SCNA_status )) 

hrd <- c4 + geom_point(aes(color = PTEN_SCNA_status),
                       position = position_jitter(w = 0.3,h = 0), alpha = 0.3) +
  geom_boxplot(aes(fill = PTEN_SCNA_status),
               outlier.size = 0.07, size = 0.5) +
  xlab("") +
  ylab("HRD levels") +
  theme_article() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 7, color = "Black"),
        #legend.position = c(0.9, 0.1),
        legend.position = "none",
        #legend.key.height = unit(2, "mm"),
        legend.title = element_text(size = 0),
        legend.box.background = element_rect(color = "white"),
        strip.text = element_text(size = 6),
        strip.text.x = element_text(angle = 0),
        strip.placement = "bottom", 
        strip.background.x=element_rect(color = NA,  fill=NA), 
        strip.background.y=element_rect(color = "black",  fill=NA)) +
  scale_fill_manual(values = color) +
  scale_color_manual(values = color) +
  stat_compare_means(comparisons = list(c("HomDel", "HemDel"), 
                                        c("HemDel","Intact"), 
                                        c("HomDel","Intact")), 
                     label.y = c(0.25, 0.30, 0.15) + 
                       stat_compare_means(label.y = 2))
hrd

d4 <- ggplot(table, mapping = aes(x = PTEN_SCNA_status, 
                                  y = Fraction_Altered, 
                                  fill = PTEN_SCNA_status )) 

gen_alt <- d4 + geom_point(aes(color = PTEN_SCNA_status),
                           position = position_jitter(w = 0.3,h = 0), alpha = 0.3) +
  geom_boxplot(aes(fill = PTEN_SCNA_status),
               outlier.size = 0.07, size = 0.5) +
  xlab("") +
  ylab("% Genome Altered") +
  theme_article() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 7, color = "Black"),
        #legend.position = c(0.9, 0.1),
        legend.position = "none",
        #legend.key.height = unit(2, "mm"),
        legend.title = element_text(size = 0),
        legend.box.background = element_rect(color = "white"),
        strip.text = element_text(size = 6),
        strip.text.x = element_text(angle = 0),
        strip.placement = "bottom", 
        strip.background.x=element_rect(color = NA,  fill=NA), 
        strip.background.y=element_rect(color = "black",  fill=NA)) +
  scale_fill_manual(values = color) +
  scale_color_manual(values = color) +
  stat_compare_means(comparisons = list(c("HomDel", "HemDel"), 
                                        c("HemDel","Intact"), 
                                        c("HomDel","Intact")), 
                     label.y = c(0.25, 0.30, 0.15) + 
                       stat_compare_means(label.y = 2))
gen_alt
ith
nsilMut
hrd

ggarrange(nsilMut, ith, hrd, gen_alt, labels = c("a", "b", "c", "d"),
          label.args = list(gp = grid::gpar(font = 30, cex =
                                              1.2)))






#reviewer 2 minor 4
fig <- ggplot(table, mapping = aes(x = PTEN_SCNA_status, 
                                  y = Intratumor_Heterogeneity, 
                                  fill = PTEN_SCNA_status )) 
intHet <- fig + geom_point(aes(color = PTEN_SCNA_status),
                       position = position_jitter(w = 0.3,h = 0), alpha = 0.3, ) +
  geom_boxplot(aes(fill = PTEN_SCNA_status),
               outlier.size = 0.07, size = 0.5) +
  xlab("") +
  ylab("ITH levels") +
  theme_article() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 7, color = "Black"),
        #legend.position = c(0.9, 0.1),
        legend.position = "none",
        #legend.key.height = unit(2, "mm"),
        legend.title = element_text(size = 0),
        legend.box.background = element_rect(color = "white"),
        strip.text = element_text(size = 6),
        strip.text.x = element_text(angle = 0),
        strip.placement = "bottom", 
        strip.background.x=element_rect(color = NA,  fill=NA), 
        strip.background.y=element_rect(color = "black",  fill=NA)) +
  scale_fill_manual(values = color) +
  scale_color_manual(values = color) +
  stat_compare_means(comparisons = list(c("Dupl", "Intact")), 
                     label.y = c(0.25, 0.30, 0.15) + 
                       stat_compare_means(label.y = 2))
intHet

figb <- ggplot(table, mapping = aes(x = PTEN_SCNA_status, 
                                   y = Aneuploidy_Score, 
                                   fill = PTEN_SCNA_status )) 
aneu <- figb + geom_point(aes(color = PTEN_SCNA_status),
                           position = position_jitter(w = 0.3,h = 0), alpha = 0.3, ) +
  geom_boxplot(aes(fill = PTEN_SCNA_status),
               outlier.size = 0.07, size = 0.5) +
  xlab("") +
  ylab("Aneuploidy Score") +
  theme_article() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 7, color = "Black"),
        #legend.position = c(0.9, 0.1),
        legend.position = "none",
        #legend.key.height = unit(2, "mm"),
        legend.title = element_text(size = 0),
        legend.box.background = element_rect(color = "white"),
        strip.text = element_text(size = 6),
        strip.text.x = element_text(angle = 0),
        strip.placement = "bottom", 
        strip.background.x=element_rect(color = NA,  fill=NA), 
        strip.background.y=element_rect(color = "black",  fill=NA)) +
  scale_fill_manual(values = color) +
  scale_color_manual(values = color) +
  stat_compare_means(comparisons = list(c("Dupl", "Intact")), 
                     label.y = c(0.25, 0.30, 0.15) + 
                       stat_compare_means(label.y = 2))
aneu

ggarrange(intHet, aneu, labels = c("a", "b"),
          label.args = list(gp = grid::gpar(font = 30, cex =
                                              1.2)))

