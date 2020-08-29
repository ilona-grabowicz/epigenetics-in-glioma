# Author: Ilona E. Grabowicz
# Date: 2019-2020

library(tidyr)
library(dplyr)
library(ggplot2)
library(data.table)
library(ggsignif)

### Load gene expression data

rna_original <- read.csv ('./input_files/RNAseq/combined_RNAs_normalised.csv', header=T, sep='\t', stringsAsFactors = F)
rna_original <- separate(rna_original,col = gene, into=c('gene', 'symbol', 'location'), sep = '_')

### Selection of genes with mean exprs>10 in all samples

gAbove10<- read.csv('GeneNamesExprAbove10.csv')

### Chromatin marks counts normalised

h3k27ac <- read.csv('./input_files/ChIPseq_ATAC_DNase/H3K27ac_coverage_quantile_norm_on_promoters_genes_selected.csv', stringsAsFactors = F, sep=',')
h3k4me3 <- read.csv('./input_files/ChIPseq_ATAC_DNase/H3K4me3_coverage_quantile_norm_on_promoters_genes_selected.csv', stringsAsFactors = F, sep=',')
atac <- read.csv('./input_files/ChIPseq_ATAC_DNase/ATACseq_coverage_quantile_norm_on_promoters_genes_selected.csv', stringsAsFactors = F, sep=',')

# Function for calculating the correlations

# what <- 'h3k4me3'
# subtype_1 <- 'PA'
# subtype_2 <- 'GB|PG'
# count_or_visualise <- 'visualise'
# read_count_thres <- 300
funkcja <- function(what, subtype_1, subtype_2, count_or_visualise, read_count_thres) { 
  my_pattern_list<-c(subtype_1, subtype_2)
  rna <- rna_original
  h3k27ac <- read.csv('./input_files/ChIPseq_ATAC_DNase/H3K27ac_coverage_quantile_norm_on_promoters_genes_selected.csv', stringsAsFactors = F, sep=',')
  h3k4me3 <- read.csv('./input_files/ChIPseq_ATAC_DNase/H3K4me3_coverage_quantile_norm_on_promoters_genes_selected.csv', stringsAsFactors = F, sep=',')
  h3k27ac$transcript <- h3k27ac$X <- h3k4me3$transcript <- h3k4me3$X <- atac$transcript <- atac$X <-NULL
  if ((subtype_1=='PA'|subtype_2=='PA') & (subtype_1=='DA'|subtype_2=='DA')){
   de <- read.csv('DeSeq2_grade1_grade23.csv', stringsAsFactors = F, header = T, sep='\t')
  }
  if ((subtype_1=='PA'|subtype_2=='PA') & (subtype_1=='GB|PG'|subtype_2=='GB|PG')){
    de <- read.csv('DeSeq2_grade1_grade4.csv', stringsAsFactors = F, header = T, sep='\t')
  }
  if ((subtype_1=='DA'|subtype_2=='DA') & (subtype_1=='GB|PG'|subtype_2=='GB|PG')){
    de <- read.csv('DeSeq2_grade23_grade4.csv', stringsAsFactors = F, header = T, sep='\t')
  }
  de <- separate(de, col = gene, into=c('gene', 'symbol', 'location'), sep = '_')
  gAboveThresh <- de$gene[de$mean_both>read_count_thres]
  de <- na.omit(de[de$FDR<0.01,])$gene
  de <- de[de%in%gAbove10$x]
  
  gene_list <- de
  gene_list <- gene_list[gene_list%in%rna$gene]
  if (what=='h3k27ac'){
    data <- h3k27ac[h3k27ac$gene%in%gene_list,]
  }
  if (what=='h3k4me3'){
    data <- h3k4me3[h3k4me3$gene%in%gene_list,]
  }
  patients <- colnames(data) [grepl(paste(my_pattern_list, collapse='|'), colnames(data))]
  rna <- rna_original[rna_original$gene%in%gene_list,c('gene', patients)]
  kor_s <- c()
  mean_expr <- c()
  genes_high_cor <- c()
  
  for (gene in gene_list){
    kor_s <- c(kor_s, cor.test(as.numeric(rna[rna$gene==gene,patients]), as.numeric(data[data$gene==gene,patients]), method='spearman')$estimate)
    mean_expr <- c(mean_expr, as.numeric(rna[rna$gene==gene,patients]))
    if (!is.na(cor.test(as.numeric(rna[rna$gene==gene,patients]), as.numeric(data[data$gene==gene,patients]), method='spearman')) & 
        cor.test(as.numeric(rna[rna$gene==gene,patients]), as.numeric(data[data$gene==gene,patients]), method='spearman')$estimate >= 0.7){
      genes_high_cor <- c(genes_high_cor, gene)
    }
  }
  
  gene_list_random <- sample(as.character(gAboveThresh), length(de)) 
  rna <- rna_original[rna_original$gene%in%gene_list_random,c('gene', patients)]
  kor_s_random <-c()
  mean_expr_random <- c()
  genes_high_cor_random <-c()
  if (what=='h3k27ac'){
    data <- h3k27ac[h3k27ac$gene%in%gene_list_random,]
  }
  if (what=='h3k4me3'){
    data <- h3k4me3[h3k4me3$gene%in%gene_list_random,]
  }
  for (gene in gene_list_random){
    kor_s_random <- c(kor_s_random, cor.test(as.numeric(rna[rna$gene==gene,patients]), as.numeric(data[data$gene==gene,patients]), method='spearman')$estimate)
    mean_expr_random <- c(mean_expr_random, as.numeric(rna[rna$gene==gene,patients]))
    if (!is.na(cor.test(as.numeric(rna[rna$gene==gene,patients]), as.numeric(data[data$gene==gene,patients]), method='spearman')) & 
        cor.test(as.numeric(rna[rna$gene==gene,patients]), as.numeric(data[data$gene==gene,patients]), method='spearman')$estimate >= 0.7){
      genes_high_cor_random <- c(genes_high_cor_random, gene)
    }
  }
  if (count_or_visualise =='count'){
    proteinatlas <- read.csv('proteinatlas.tsv', sep = '\t', stringsAsFactors = F)
    proteinatlas_glioma <- proteinatlas$Ensembl[grepl('Glioma', proteinatlas$Prognostic.p.value)]
    print(length(genes_high_cor))
    genes_corr_prognostic_ENSEMBL <-c(genes_high_cor[genes_high_cor %in% proteinatlas_glioma])
    genes_corr_prognostic_Symbol <- rna_original$symbol[rna_original$gene %in% genes_corr_prognostic_ENSEMBL]
    print(genes_corr_prognostic_Symbol)
    print(sum(genes_high_cor %in% proteinatlas_glioma) / length(de) *100) # 1.28)
    print(sum(genes_high_cor_random %in% proteinatlas_glioma) / length(de) *100) # 0.21
  }
  else{return(list('kor_s'=kor_s, 'kor_s_random'=kor_s_random, 'mean_expr'=mean_expr, 'mean_expr_random' = mean_expr_random))}
}

# To calculate the mean of random shares run this e.g. 10x
# 'h3k4me3', 'PA', 'DA': 
cor_de_1_23_h3k4me3<-funkcja('h3k4me3', 'PA', 'DA', 'count', 300) #  1.28 
(0.42+0.47+0.25+0.34+0.42+0.42+0.38+0.64+0.34+0.29)/10 # 0.39

# 'h3k27ac', 'PA', 'DA': 
cor_de_1_23_h3k27ac<-funkcja('h3k27ac', 'PA', 'DA', 'count', 300) # 1.32
(0.51+0.38+0.47+0.68+0.42+0.55+0.47+0.51+0.42+0.51)/10 # 0.49

# 'h3k4me3', 'DA', 'GB|PG': 
cor_de_23_4_h3k4me3<-funkcja('h3k4me3', 'DA', 'GB|PG', 'count', 100) # 3.06
(0.00+0.00+0.00+0.00+1.02+0.00+0.00+0.00+0.00+0.00)/10 # 0.102

# 'h3k27ac', 'DA', 'GB|PG': 
cor_de_23_4_h3k27ac<-funkcja('h3k27ac', 'DA', 'GB|PG', 'count', 100) # 5.1
(1.02+1.02+0.00+0.00+1.02+1.02+2.04+1.02+0.00+0.00)/10 # 0.7

# 'h3k4me3', 'PA', 'GB|PG': 
cor_de_1_4_h3k4me3<-funkcja('h3k4me3', 'PA', 'GB|PG', 'count',200) # 0.24
(0.12+0.12+0.00+0.09+0.12+0.06+0.12+0.18+0.12+0.06)/10 # 0.09

# 'h3k27ac', 'PA', 'GB|PG': 
cor_de_1_4_h3k27ac<-funkcja('h3k27ac', 'PA', 'GB|PG', 'count',200) # 0.60
(0.33+0.30+0.36+0.18+0.24+0.36+0.30+0.33+0.33+0.36)/10 # 0.30

wilcox.test(c(1.28, 1.32, 3.06, 5.1, 0.24, 0.6), c(0.39, 0.49, 0.10, 0.7, 0.09, 0.3)) # p-value = 0.04
mean(1.28, 1.32, 3.06, 5.1, 0.24, 0.6) #1.28
sd(c(1.28, 1.32, 3.06, 5.1, 0.24, 0.6)) #1.83
mean(c(0.39, 0.49, 0.10, 0.7, 0.09, 0.3)) #0.34
sd(c(0.39, 0.49, 0.10, 0.7, 0.09, 0.3)) #0.23

######### Visualisations:

# PA vs DA
cor_de_1_23_h3k4me3<-funkcja('h3k4me3', 'PA', 'DA', 'visualise', 300)
cor_de_1_23_h3k27ac<-funkcja('h3k27ac', 'PA', 'DA', 'visualise', 300)

mean(cor_de_1_23_h3k4me3[[3]]) 
mean(cor_de_1_23_h3k4me3[[4]])  
mean(cor_de_1_23_h3k27ac[[3]]) 
mean(cor_de_1_23_h3k27ac[[4]])  
median(cor_de_1_23_h3k4me3[[1]]) #0.40
median(cor_de_1_23_h3k4me3[[2]]) #0.25
median(cor_de_1_23_h3k27ac[[1]]) #0.57
median(cor_de_1_23_h3k27ac[[2]]) #0.33

df <- data.frame("korelacje" = c(cor_de_1_23_h3k4me3[[1]], cor_de_1_23_h3k27ac[[1]], cor_de_1_23_h3k4me3[[2]], cor_de_1_23_h3k27ac[[2]]), 
                 "variable" = c(rep(paste0("H3K4me3\n(n=", length(cor_de_1_23_h3k4me3[[1]]), ")"), length(cor_de_1_23_h3k4me3[[1]])),
                                rep(paste0("H3K27ac\n(n=", length(cor_de_1_23_h3k27ac[[1]]), ")"), length(cor_de_1_23_h3k27ac[[1]])),
                                rep(paste0("random genes\n H3K4me3\n(n=", length(cor_de_1_23_h3k4me3[[2]]), ")"), length(cor_de_1_23_h3k4me3[[2]])),
                                rep(paste0("random genes\n H3K27ac\n(n=", length(cor_de_1_23_h3k27ac[[2]]), ")"), length(cor_de_1_23_h3k27ac[[2]]))
                 )
)

df$variable <- factor(df$variable,
                      levels = c(paste0("H3K4me3\n(n=", length(cor_de_1_23_h3k4me3[[1]]),")"),
                                 paste0("random genes\n H3K4me3\n(n=", length(cor_de_1_23_h3k4me3[[2]]),")"),
                                 paste0("H3K27ac\n(n=", length(cor_de_1_23_h3k27ac[[1]]), ")"),
                                 paste0("random genes\n H3K27ac\n(n=", length(cor_de_1_23_h3k27ac[[2]]),")")
                                 ))#,ordered = TRUE)

ggplot(df, aes(x=variable, y=korelacje)) +
  geom_boxplot(color=c('coral', 'coral', 'dodgerblue', 'dodgerblue'),
               fill=c('white', 'grey', 'white', 'grey')) +
  xlab('') +
  scale_y_continuous('Spearman correlation rho', breaks = seq(-1, 1, by=0.2), limits = c(-1, 1))+
  theme(plot.title = element_text(size = 14)) +
  theme(text = element_text(size=15)) +
  geom_signif(annotation='***', xmin=1, xmax=2, y_position=1) +
  geom_signif(annotation='***', xmin=3, xmax=4, y_position=1)
ggsave('DE_1_23_both_epi.png', width = 6.5, height = 5)

wilcox.test(cor_de_1_23_h3k4me3[[3]], cor_de_1_23_h3k4me3[[4]])$p.value # 0
wilcox.test(cor_de_1_23_h3k27ac[[3]], cor_de_1_23_h3k27ac[[4]])$p.value # 0


# PA vs GB|PG
cor_de_1_4_h3k4me3<-funkcja('h3k4me3', 'PA', 'GB|PG', 'visualise', 100)
cor_de_1_4_h3k27ac<-funkcja('h3k27ac', 'PA', 'GB|PG', 'visualise', 100)

mean(cor_de_1_4_h3k4me3[[3]]) 
mean(cor_de_1_4_h3k4me3[[4]])  
mean(cor_de_1_4_h3k27ac[[3]]) 
mean(cor_de_1_4_h3k27ac[[4]]) 
median(cor_de_1_4_h3k4me3[[1]]) #0.34
median(cor_de_1_4_h3k4me3[[2]]) #0.22
median(cor_de_1_4_h3k27ac[[1]]) #0.49
median(cor_de_1_4_h3k27ac[[2]]) #0.32

df <- data.frame("korelacje" = c(cor_de_1_4_h3k4me3[[1]], cor_de_1_4_h3k27ac[[1]], cor_de_1_4_h3k4me3[[2]], cor_de_1_4_h3k27ac[[2]]), 
                 "variable" = c(rep(paste0("H3K4me3\n(n=", length(cor_de_1_4_h3k4me3[[1]]), ")"), length(cor_de_1_4_h3k4me3[[1]])),
                                rep(paste0("H3K27ac\n(n=", length(cor_de_1_4_h3k27ac[[1]]), ")"), length(cor_de_1_4_h3k27ac[[1]])),
                                rep(paste0("random genes\n H3K4me3\n(n=", length(cor_de_1_4_h3k4me3[[2]]), ")"), length(cor_de_1_4_h3k4me3[[2]])),
                                rep(paste0("random genes\n H3K27ac\n(n=", length(cor_de_1_4_h3k27ac[[2]]), ")"), length(cor_de_1_4_h3k27ac[[2]]))
                 )
)

df$variable <- factor(df$variable,
                      levels = c(paste0("H3K4me3\n(n=", length(cor_de_1_4_h3k4me3[[1]]),")"),
                                 paste0("random genes\n H3K4me3\n(n=", length(cor_de_1_4_h3k4me3[[2]]),")"),
                                 paste0("H3K27ac\n(n=", length(cor_de_1_4_h3k27ac[[1]]), ")"),
                                 paste0("random genes\n H3K27ac\n(n=", length(cor_de_1_4_h3k27ac[[2]]),")")
                      ))#,ordered = TRUE)

ggplot(df, aes(x=variable, y=korelacje)) +
  geom_boxplot(color=c('coral', 'coral', 'dodgerblue', 'dodgerblue'),
               fill=c('white', 'grey', 'white', 'grey')) +
  xlab('') +# Gene sets ylab() +
  scale_y_continuous('Spearman correlation rho', breaks = seq(-1, 1, by=0.2), limits = c(-1, 1))+
  theme(plot.title = element_text(size = 14)) +
  theme(text = element_text(size=15)) +
  geom_signif(annotation='***', xmin=1, xmax=2, y_position=1) +
  geom_signif(annotation='***', xmin=3, xmax=4, y_position=1)
ggsave('DE_1_4_both_epi.png', width = 6.5, height = 5)


wilcox.test(cor_de_1_4_h3k4me3[[3]], cor_de_1_4_h3k4me3[[4]])$p.value # 0
wilcox.test(cor_de_1_4_h3k27ac[[3]], cor_de_1_4_h3k27ac[[4]])$p.value # 0

# DA vs GB|PG
cor_de_23_4_h3k4me3<-funkcja('h3k4me3', 'DA', 'GB|PG', 'visualise', 200)
cor_de_23_4_h3k27ac<-funkcja('h3k27ac', 'DA', 'GB|PG', 'visualise', 200)

mean(cor_de_23_4_h3k4me3[[3]]) 
mean(cor_de_23_4_h3k4me3[[4]])  
mean(cor_de_23_4_h3k27ac[[3]]) 
mean(cor_de_23_4_h3k27ac[[4]]) 
median(cor_de_23_4_h3k4me3[[1]]) #0.42
median(cor_de_23_4_h3k4me3[[2]]) #0.23
median(cor_de_23_4_h3k27ac[[1]]) #0.66
median(cor_de_23_4_h3k27ac[[2]]) #0.42

df <- data.frame("korelacje" = c(cor_de_23_4_h3k4me3[[1]], cor_de_23_4_h3k27ac[[1]], cor_de_23_4_h3k4me3[[2]], cor_de_23_4_h3k27ac[[2]]), 
                 "variable" = c(rep(paste0("H3K4me3\n(n=", length(cor_de_23_4_h3k4me3[[1]]), ")"), length(cor_de_23_4_h3k4me3[[1]])),
                                rep(paste0("H3K27ac\n(n=", length(cor_de_23_4_h3k27ac[[1]]), ")"), length(cor_de_23_4_h3k27ac[[1]])),
                                rep(paste0("random genes\n H3K4me3\n(n=", length(cor_de_23_4_h3k4me3[[2]]), ")"), length(cor_de_23_4_h3k4me3[[2]])),
                                rep(paste0("random genes\n H3K27ac\n(n=", length(cor_de_23_4_h3k27ac[[2]]), ")"), length(cor_de_23_4_h3k27ac[[2]]))
                 )
)

df$variable <- factor(df$variable,
                      levels = c(paste0("H3K4me3\n(n=", length(cor_de_23_4_h3k4me3[[1]]),")"),
                                 paste0("random genes\n H3K4me3\n(n=", length(cor_de_23_4_h3k4me3[[2]]),")"),
                                 paste0("H3K27ac\n(n=", length(cor_de_23_4_h3k27ac[[1]]), ")"),
                                 paste0("random genes\n H3K27ac\n(n=", length(cor_de_23_4_h3k27ac[[2]]),")")
                      ))#,ordered = TRUE)

ggplot(df, aes(x=variable, y=korelacje)) +
  geom_boxplot(color=c('coral', 'coral', 'dodgerblue', 'dodgerblue'),
               fill=c('white', 'grey', 'white', 'grey')) +
  xlab('') +
  scale_y_continuous('Spearman correlation rho', breaks = seq(-1, 1, by=0.2), limits = c(-1, 1))+
  theme(plot.title = element_text(size = 14)) +
  theme(text = element_text(size=15)) +
  geom_signif(annotation='***', xmin=1, xmax=2, y_position=1) +
  geom_signif(annotation='***', xmin=3, xmax=4, y_position=1)
ggsave('DE_23_4_both_epi.png', width = 6.5, height = 5)

wilcox.test(cor_de_23_4_h3k4me3[[3]], cor_de_23_4_h3k4me3[[4]])$p.value # 2.17155e-246
wilcox.test(cor_de_23_4_h3k27ac[[3]], cor_de_23_4_h3k27ac[[4]])$p.value # 9.147805e-187


# Differences in correlations between real DE and random genes
# H3K4me3:
mean(0.40-0.25, 0.34-0.22, 0.42-0.23) #0.15
# H3K27ac:
mean(0.57-0.33, 0.49-0.32, 0.66-0.42) #0.24
