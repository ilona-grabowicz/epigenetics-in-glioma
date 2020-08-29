# Author: Ilona E. Grabowicz
# Date: 2019-2020

library(enrichR)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DOSE)
library(GO.db)
library(clusterProfiler)
library(org.Hs.eg.db)

########## Finding enriched enriched GO terms for different sets of genes, such as sets from enriched glioma TADs

enriched_tads <- read.csv('./TADs/counting_DE_TADS.csv')
enriched_tads <- enriched_tads$Var1[enriched_tads$Freq>=3]
tads_with_genes <- read.csv('./TADs/TADs_with_genes_ENSEMBL_and_FC_RNASeq_grade1_grade23.csv')
#tad <- 2007
threshold <- 0.01
for (tad in enriched_tads){
  genes <- tads_with_genes[tads_with_genes$TAD_nr==tad, grepl('gene', colnames(tads_with_genes))]
  genes <- genes[!is.na(genes)]
  input_genes <- genes
  
  ### EnrichR functional annotation
  
  gene.df <- bitr(input_genes, fromType = "ENSEMBL", #names(geneList_FC)
                  toType = c("ENTREZID", "SYMBOL"),
                  OrgDb = org.Hs.eg.db)
  databases <- c("GO_Biological_Process_2018")#, "GO_Molecular_Function_2018", "TRANSFAC_and_JASPAR_PWMs", "KEGG_2019_Human", "Panther_2016",
                # "InterPro_Domains_2019")
  #db <- "GO_Biological_Process_2018"
  for (db in databases) {
    enrichr_result<-enrichr(gene.df$SYMBOL, databases = db)
    write.csv(enrichr_result[[1]], file = paste0('./TADs/gene_enrichments/enrichr_TAD_',tad,'_', db, '.csv'))
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
      ggtitle(paste0('TAD ', tad)) +
      labs(size="Count of genes\ninvolved\nin the term", color='p-value') +
      scale_color_gradient(low="red", high="blue")  +
      theme(text = element_text(size=15)) +
      guides(size = guide_legend(order = 1), 
              colour = guide_legend(order = 2))
    ggsave(paste0('./TADs/gene_enrichments/TAD_', tad, '_', db, '_', threshold, '.png'), width = 12, height = 5)
  }
}

# Enrichments for other sets of genes:
# input_genes <- c('SNX10', 'PHGDH', 'QPRT', 'CPQ', 'AEBP1', 'ADAMTSL4', 'UBXN10', 'EN2', 'CARD19', 'HTRA1', 'ISG20', 'GPC5', 'MKRN3', 'ACHE', 
#                  'MXRA5', 'KCNN4', 'RARRES2', 'LOXL2', 'MSTN', 'ETNK2', 'RNF175', 'FBXO27', 'CYGB', 'RPL39L', 'ETNPPL', 'FOSL1', 'HIST3H2A', 'CASC10', 
#                  'EMP2', 'FBXO17')
# input_genes <- c('CPQ', 'UBXN10', 'MXRA5', 'RARRES2', 'SYT5', 'CYGB', 'CEND1', 'EMP2')
# input_genes <- c('ISG20', 'GLUD1', 'PDGFA')
# input_genes <- c('MTHFD2', 'SNX10', 'PHGDH', 'QPRT', 'CPQ', 'AEBP1', 'LOXL1', 'DNAJA4', 'CTF1', 'CRELD1', 'EN2', 'CARD19', 'SMPD1', 'ISG20', 'MKRN3', 
#                  'ACHE', 'RARRES2', 'LOXL4', 'MSTN', 'ETNK2', 'RNF175', 'GLUD1', 'FBXO27', 'CYGB', 'ETNPPL', 'C11orf24', 'GRIK1', 'HIST3H2A', 
#                  'RP11-195F19.5', 'CASC10', 'FBXO17')
# input_genes <- c('RWDD2A', 'SPAG4', 'QPRT', 'PDLIM2', 'HS3ST2', 'REEP2', 'ITPKA', 'RABGEF1', 'UBXN10', 'HTRA1', 'ACHE', 'IMPDH1', 'RARRES2', 'STAC2', 'THY1', 'STC1', 'PRSS12',
#                  'C11orf24', 'CEND1', 'EMP2')
# input_genes <- c('ISG20', 'GLUD1', 'STC1', 'FOSL1', 'PDGFA')

### EnrichR functional annotation
threshold <- 0.005
draw <- function(input) {
  gene.df <- bitr(input, fromType = "SYMBOL", #names(geneList_FC)
                  toType = c("ENTREZID", "SYMBOL"),
                  OrgDb = org.Hs.eg.db)
  #databases <- c("GO_Biological_Process_2018")#, "GO_Molecular_Function_2018", "TRANSFAC_and_JASPAR_PWMs", "KEGG_2019_Human", "Panther_2016",
  #  "InterPro_Domains_2019")
  db <- "GO_Biological_Process_2018"
  #for (db in databases) {
  enrichr_result<-enrichr(gene.df$SYMBOL, databases = db)
  write.csv(enrichr_result[[1]], file = paste0('./TFs/enrichr_db.csv'))
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
    labs(size="Count of genes\ninvolved\nin the term", color='p-value') +
    scale_color_gradient(low="red", high="blue") +
    theme(text = element_text(size=15)) 
  
  #print(g)
  # return(g)
  #  }
}

# Plots for glioma TADs:

input_genes <- c('SNX10', 'PHGDH', 'QPRT', 'CPQ', 'AEBP1', 'ADAMTSL4', 'UBXN10', 'EN2', 'CARD19', 'HTRA1', 'ISG20', 'GPC5', 'MKRN3', 'ACHE', 
                 'MXRA5', 'KCNN4', 'RARRES2', 'LOXL2', 'MSTN', 'ETNK2', 'RNF175', 'FBXO27', 'CYGB', 'RPL39L', 'ETNPPL', 'FOSL1', 'HIST3H2A', 'CASC10', 
                 'EMP2', 'FBXO17')
threshold <- 0.005
g <- draw(input_genes)
if (length(g$data$Term) > 15) {
  threshold <- 0.0025
  g <- draw(input_genes)
}
g <- g + ggtitle(paste0('H3K4me3, PA vs DA'))
print(g)
ggsave(paste0('./TFs/GO_BP_H3K4me3_PA_DA_', threshold,'.png'), width = 12, height = 5)



input_genes <- c('CPQ', 'UBXN10', 'MXRA5', 'RARRES2', 'SYT5', 'CYGB', 'CEND1', 'EMP2')
threshold <- 0.005
g <- draw(input_genes)
if (length(g$data$Term) > 17) {
  threshold <- 0.0025
  g <- draw(input_genes)
}
g <- g + ggtitle(paste0('H3K4me3, PA vs GB/PG'))
print(g)
ggsave(paste0('./TFs/GO_BP_H3K4me3_PA_GB_PG_', threshold,'.png'), width = 12, height = 5)



input_genes <- c('ISG20', 'GLUD1', 'PDGFA')
threshold <- 0.005
g <- draw(input_genes)
if (length(g$data$Term) > 17) {
  threshold <- 0.0025
  g <- draw(input_genes)
}
g <- g + ggtitle(paste0('H3K4me3, DA vs GB/PG'))
print(g)
ggsave(paste0('./TFs/GO_BP_H3K4me3_DA_GB_PG_', threshold,'.png'), width = 12, height = 5)



input_genes <- c('MTHFD2', 'SNX10', 'PHGDH', 'QPRT', 'CPQ', 'AEBP1', 'LOXL1', 'DNAJA4', 'CTF1', 'CRELD1', 'EN2', 'CARD19', 'SMPD1', 'ISG20', 'MKRN3', 
                 'ACHE', 'RARRES2', 'LOXL4', 'MSTN', 'ETNK2', 'RNF175', 'GLUD1', 'FBXO27', 'CYGB', 'ETNPPL', 'C11orf24', 'GRIK1', 'HIST3H2A', 
                 'RP11-195F19.5', 'CASC10', 'FBXO17')
threshold <- 0.005
g <- draw(input_genes)
if (length(g$data$Term) > 15) {
  threshold <- 0.0025
  g <- draw(input_genes)
}
g <- g + ggtitle(paste0('H3K27ac, PA vs DA'))
print(g)
ggsave(paste0('./TFs/GO_BP_H3K27ac_PA_DA_', threshold,'.png'), width = 12, height = 5)




input_genes <- c('RWDD2A', 'SPAG4', 'QPRT', 'PDLIM2', 'HS3ST2', 'REEP2', 'ITPKA', 'RABGEF1', 'UBXN10', 'HTRA1', 'ACHE', 'IMPDH1', 'RARRES2', 'STAC2', 'THY1', 'STC1', 'PRSS12',
                 'C11orf24', 'CEND1', 'EMP2')
threshold <- 0.005
g <- draw(input_genes)
if (length(g$data$Term) > 15) {
  threshold <- 0.0025
  g <- draw(input_genes)
}
g <- g + ggtitle(paste0('H3K27ac, PA vs GB/PG'))
print(g)
ggsave(paste0('./TFs/GO_BP_H3K27ac_PA_GB_PG_', threshold,'.png'), width = 12, height = 5)



input_genes <- c('ISG20', 'GLUD1', 'STC1', 'FOSL1', 'PDGFA')
threshold <- 0.005
g <- draw(input_genes)
if (length(g$data$Term) > 17) {
  threshold <- 0.0025
  g <- draw(input_genes)
}
g <- g + ggtitle(paste0('H3K27ac, DA vs GB/PG'))
print(g)
ggsave(paste0('./TFs/GO_BP_H3K27ac_DA_GB_PG_', threshold,'.png'), width = 12, height = 5)

