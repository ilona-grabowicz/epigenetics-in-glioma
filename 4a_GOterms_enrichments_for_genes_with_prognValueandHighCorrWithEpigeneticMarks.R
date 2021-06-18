library(tidyr)
library(dplyr)
library(ggplot2)
library(data.table)
library(enrichR)
library(DOSE)
library(GO.db)
library(clusterProfiler)
library(org.Hs.eg.db)

#setwd('/mnt/chr3/People/Ilona/cancer/2019_10_29_code_for_Github/')

######### Load gene expression data

rna_original <- read.csv ('./input_files/RNAseq/combined_RNAs_normalised.csv', header=T, sep='\t', stringsAsFactors = F)
rna_original <- separate(rna_original,col = gene, into=c('gene', 'symbol', 'location'), sep = '_')

### All genes with mean exprs>10 in all samples

gAbove10<- read.csv('GeneNamesExprAbove10.csv')

### Chromatin marks counts normalised

h3k27ac <- read.csv('./input_files/ChIPseq_ATAC_DNase/H3K27ac_coverage_quantile_norm_on_promoters_genes_selected.csv', stringsAsFactors = F, sep=',')
h3k4me3 <- read.csv('./input_files/ChIPseq_ATAC_DNase/H3K4me3_coverage_quantile_norm_on_promoters_genes_selected.csv', stringsAsFactors = F, sep=',')

######### Filtering out the genes:

what <- 'h3k27ac'
subtype_1 <- 'PA'
subtype_2 <- 'GB|PG'
read_count_thres <- 300

my_pattern_list<-c(subtype_1, subtype_2)
rna <- rna_original
h3k27ac$transcript <- h3k27ac$X <- h3k4me3$transcript <- h3k4me3$X <-NULL
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

# Filtering out only genes with active expression levels:
gAboveThresh <- de$gene[de$mean_both>read_count_thres]

# Filtering out only genes significantly differential:
de <- na.omit(de[de$FDR<0.01,])$gene
de <- de[de%in%gAbove10$x]
gene_list <- de

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
  # Filtering out only genes with high correlation with epigenetic mark:
  if (!is.na(cor.test(as.numeric(rna[rna$gene==gene,patients]), as.numeric(data[data$gene==gene,patients]), method='spearman')) & 
      cor.test(as.numeric(rna[rna$gene==gene,patients]), as.numeric(data[data$gene==gene,patients]), method='spearman')$estimate >= 0.7){
    genes_high_cor <- c(genes_high_cor, gene)
  }
}
# Filtering out only genes with prognostic value for glioma:
proteinatlas <- read.csv('proteinatlas.tsv', sep = '\t', stringsAsFactors = F)
proteinatlas_glioma <- proteinatlas$Ensembl[grepl('Glioma', proteinatlas$Prognostic.p.value)]
genes_corr_prognostic_ENSEMBL <-c(genes_high_cor[genes_high_cor %in% proteinatlas_glioma])
genes_corr_prognostic_Symbol <- rna_original$symbol[rna_original$gene %in% genes_corr_prognostic_ENSEMBL]
# End list of genes:
print(genes_corr_prognostic_Symbol)


#########  Plotting the enriched GO terms:
input_genes <- genes_corr_prognostic_Symbol
threshold <- 0.005
draw <- function(input) {
  gene.df <- bitr(input, fromType = "SYMBOL", #names(geneList_FC)
                  toType = c("ENTREZID", "SYMBOL"),
                  OrgDb = org.Hs.eg.db)
  db <- "GO_Biological_Process_2018"
  enrichr_result<-enrichr(gene.df$SYMBOL, databases = db)
  data_to_plot <- data.frame(enrichr_result[[1]])
  # PLOT:
  data_to_plot <- separate(data = data_to_plot, col=Overlap, into=c('nr_genes', 'nr_genes_in_term'), sep='/')
  data_to_plot$Percentage <- as.numeric(data_to_plot$nr_genes)/as.numeric(data_to_plot$nr_genes_in_term)*100
  data_to_plot$nr_genes <- as.integer(data_to_plot$nr_genes)
  data_to_plot$Term <- factor(as.character(data_to_plot$Term), levels = data_to_plot$Term[order(data_to_plot$Percentage)])
  data_to_plot <- data_to_plot[data_to_plot$P.value<threshold,]
  g <- ggplot(data_to_plot, aes(x=Percentage, y=Term, color=P.value, size = factor(nr_genes))) +
    geom_point() +
    xlab('Percentage of DE genes involved in the term') +
    ylab( db) +
    #ggtitle(paste0('H3K27ac, PA vs GB/PG')) +
    labs(size="Count of genes\ninvolved\nin the term", color='p-value') +
    scale_color_gradient(low="red", high="blue") +
    theme(text = element_text(size=14)) 
}
g <- draw(input_genes)
g <- g + ggtitle(paste0('H3K27ac, PA vs GBM/pGBM'))
print(g)


# Saving Fig 1D:
ggsave(paste0('./GO_BP_H3K27ac_PA_GBMpGBM_', threshold,'.svg'), width = 12, height = 5)
write.csv(x = data_to_plot, file = 'GO_terms_H3K27ac_PAvsGBMpGBM.csv')
