# Author: Ilona E. Grabowicz
# Date: 2019-2020

library(ggplot2)
library(reshape)
library(dplyr)
library(tidyr)

############## Finding which genes have bivalent chromatin in grade IV samples (GB and PG):

h3k27me3 <- read.csv('./input_files/ChIPseq_ATAC_DNase/H3K27me3_coverage_on_promoters.tsv', stringsAsFactors = F, sep='\t', header=T)
h3k27me3$X<-NULL
h3k4me3 <- read.csv('./input_files/ChIPseq_ATAC_DNase/H3K4me3_coverage_on_promoters.tsv', stringsAsFactors = F, sep='\t', header=T)
h3k4me3$X <- NULL

grade4 <- colnames(h3k4me3)[as.logical(grepl('GB', colnames(h3k4me3)) + grepl('PG', colnames(h3k4me3)))]
grade4 <- grade4[grade4 %in% colnames(h3k27me3)]

grade4_df <- data.frame("gene_ID" = h3k27me3$gene_ID)
grade4_df$h3k4me3_mean <- apply(h3k4me3[,colnames(h3k4me3)%in%grade4], 1, mean) 
grade4_df$h3k27me3_mean <- apply(h3k27me3[,colnames(h3k27me3)%in%grade4], 1, mean) 

tmp<- arrange(grade4_df, desc(grade4_df$h3k4me3_mean))
k4_top<-as.character(tmp$gene_ID[1:1000])
tmp<- arrange(grade4_df, desc(grade4_df$h3k27me3_mean))
k27_top<-as.character(tmp$gene_ID[1:1000])
top1000<- intersect(k4_top, k27_top)
write.csv(top1000, 'Suppl_Table_5_54_DEGs_with_bivalent_chromatin_in_GB.csv')

deseq<- read.csv('./DeSeq2_grade1_grade4.csv', stringsAsFactors = F, sep='\t', header=T)
deseq <- separate(deseq,col = gene, into=c('ENSEMBL', 'symbol', 'location'), sep = '_') 
deseq <- separate(deseq, col=location, into=c('chromosome', 'rest'), sep=':')
deseq <- separate(deseq, col=rest, into=c('start', 'stop'))
deseq <- separate(deseq, col=chromosome, into=c('chromosome', 'x'))
deseq$x<-NULL

# Narrowing down deseq result to top 1000 H3K27me3 and H3K4me3 genes
deseq<- deseq[deseq$ENSEMBL%in%top1000,]
deseq$tad<- NA
TADs <- read.csv('./input_files/TADborders-Fetal_brain_D_cortical_plates_Won_et_al_2016.bed', stringsAsFactors = F)
list_of_genes <-list()
list_of_TADs <-c()
 
# Assigning TADs to DE genes:
for (i in 1:nrow(deseq)){ #nrow(deseq)
  for(j in 1:nrow(TADs)){
    
    if (deseq$chromosome[i] == TADs$chr[j]) {
      if (as.numeric(deseq$start[i])>= as.numeric(TADs$start[j]) & as.numeric(deseq$stop[i])<= as.numeric(TADs$stop[j])){
        deseq$tad[i] <- TADs$TAD_nr[j]
        list_of_genes[[deseq$ENSEMBL[i]]]<- TADs$TAD_nr[j]
        list_of_TADs <- c(list_of_TADs, TADs$TAD_nr[j])
        cat(deseq$ENSEMBL[i], TADs$TAD_nr[j], '\n')
      } 
    }
  }
}

regulated_TADs_file <- read.csv ('./TADs/counting_DE_TADS.csv', stringsAsFactors = F)
regulated_TADs <- regulated_TADs_file$Var1[regulated_TADs_file$Freq>=2]
regulated_TADs <- regulated_TADs_file$Var1

intersect(list_of_TADs, regulated_TADs)
phyper(length(intersect(list_of_TADs, regulated_TADs)), length(regulated_TADs) , (3165 - length(regulated_TADs)) , length(top1000) , lower.tail = F)

cat('gene\tTAD_nr\n', file = 'genes_in_glioma_TADs.csv')
for (el in names(list_of_genes)){
  cat(paste0(el, '\t', list_of_genes[[el]], '\n'), file = 'genes_in_glioma_TADs.csv', append = T)
} 

# Plot of functions of DEGs with bivalent chromatin
library(enrichR)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DOSE)
library(GO.db)
library(clusterProfiler)
library(org.Hs.eg.db)

input_genes <- top1000
threshold <- 0.001
### EnrichR functional annotation

gene.df <- bitr(input_genes, fromType = "ENSEMBL", #names(geneList_FC)
                toType = c("ENTREZID", "SYMBOL"),
                OrgDb = org.Hs.eg.db)
# databases <- c("GO_Biological_Process_2018", "GO_Molecular_Function_2018", "TRANSFAC_and_JASPAR_PWMs", "KEGG_2019_Human", "Panther_2016",
#                "InterPro_Domains_2019")
db <- "GO_Biological_Process_2018"
#for (db in databases) {
enrichr_result<-enrichr(gene.df$SYMBOL, databases = db)

data_to_plot <- data.frame(enrichr_result[[1]])
# PLOT:
data_to_plot <- separate(data = data_to_plot, col=Overlap, into=c('nr_genes', 'nr_genes_in_term'), sep='/')
data_to_plot$Percentage <- as.numeric(data_to_plot$nr_genes)/as.numeric(data_to_plot$nr_genes_in_term)*100
data_to_plot$nr_genes <- as.integer(data_to_plot$nr_genes)
data_to_plot$Term <- factor(as.character(data_to_plot$Term), levels = data_to_plot$Term[order(data_to_plot$Percentage)])
data_to_plot <- data_to_plot[data_to_plot$P.value<threshold,]

ggplot(data_to_plot, aes(x=Percentage, y=Term, color=P.value, size = factor(nr_genes))) +
  geom_point() +
  xlab('Percentage of DE genes involved in the term') +
  ylab('') +
  #ggtitle(paste0('TAD ', tad)) +
  labs(size="Count of genes\ninvolved\nin the term", color='p-value') +
  scale_color_gradient(low="red", high="blue")  +
  theme(text = element_text(size=15)) +
  guides(size = guide_legend(order = 1), 
         colour = guide_legend(order = 2))

ggsave(paste0('./TADs/gene_enrichments/genes_wth_bivalent_chromatin.png'), height =4, width=14)
#}