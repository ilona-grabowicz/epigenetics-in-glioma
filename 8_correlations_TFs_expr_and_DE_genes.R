# Author: Ilona E. Grabowicz
# Date: 2019-2020

library(tidyr)
library(dplyr)
library(data.table)
library(ggplot2)
library(ggrepel)
library(corrplot)
library(magrittr)

################## Finding which TFs have highest correlations of expressions with DEGs:

comparison <- 'grade1_grade4'

DESeq <- read.csv(paste0('DeSeq2_', comparison, '.csv'), sep='\t')
DESeq <- separate(DESeq,col = gene, into=c('ENSEMBL', 'symbol', 'location'), sep = '_') 
DESeq<- na.omit(DESeq[DESeq$FDR<0.05,])

exprs <- read.csv('./input_files/RNAseq/combined_RNAs_normalised.csv', sep='\t', header=T, stringsAsFactors = F)
exprs <- separate(exprs,col = gene, into=c('ENSEMBL', 'symbol', 'location'), sep = '_') 
regulated_TADs <- read.csv('./TADs/counting_DE_TADS.csv')
regulated_TADs <- regulated_TADs$Var1[regulated_TADs$Freq>=3]

repeated_enriched_TFs <- setNames(data.frame(matrix(ncol = 3, nrow = 0), stringsAsFactors=F), c("enriched_TF_above_random", "TAD_nr", "median_correlation"))
# repeated_enriched_TFs jest do tego zeby zobaczyc ktore TFs maja w roznych TADach korelacje wyzsza niz randomowa
mapping <- read.csv('./input_files/mapping_TFs_not_present.csv', stringsAsFactors = F)

for (reg_tad in regulated_TADs) {
  reg_tad <- as.character(reg_tad)
  print(reg_tad)
  TFs_data <- read.csv(paste0('./TADs/enriched_TFs/TAD',reg_tad,'.txt'), sep='\t', header=T, stringsAsFactors = F) # These files are prepared by hand from David website
  if (nrow(TFs_data)>0){
    TFs_data<-TFs_data[TFs_data$FDR<=0.01,]
    if (nrow(TFs_data)>0){
      TFs_enriched <- TFs_data$Term
      #tf <- 'NRSF'
      koreluj <- function(tf) {
        TF_David <- tf
        TF_exprs <- mapping$exprs_name[mapping$David_name==tf]
        enriched_genes <- TFs_data[TFs_data$Term==TF_David,'Genes']
        enriched_genes <- strsplit(enriched_genes, ', ')
        z<-data.frame('deseq' = enriched_genes[[1]])
        enriched_genes2<- as.character(z$deseq[z$deseq %in% DESeq$ENSEMBL])
        x <- as.numeric(exprs[exprs$symbol==TF_exprs,c(4:ncol(exprs))])
        cors <- c()
        for (gene in enriched_genes2) {
          if (gene %in% exprs$ENSEMBL){
            y <- as.numeric(exprs[exprs$ENSEMBL==gene,c(4:ncol(exprs))])
            cors<- c(cors, cor.test(x,y, method='spearman')$estimate)}
        }
        return(median(cors))
      }
      
      result <- setNames(data.frame(matrix(ncol = 3, nrow = 0), stringsAsFactors=F), c("TF_David", "TF_exprs", "median_correlation"))
      
      for (tf in TFs_enriched) {
        if (tf %in% mapping$David_name){
          #print(tf)
          cors <- koreluj(tf)
          result <- rbind(result, data.frame('TF_David'=tf, 'TF_exprs'=mapping$exprs_name[mapping$David_name==tf], 'median_correlation'=cors))
        }
      }
      
      write.csv(arrange(result, desc(median_correlation)), paste0('./TADs/enriched_TFs/', reg_tad, '_TFs_correlations_onlyDE.csv'))
      
      
      ### PERMUTATIONS
      #tf <- 'NRSF'
      rysuj_random <- function(tf) {
        TF_David <- tf
        TF_exprs <- mapping$exprs_name[mapping$David_name==tf]
        enriched_genes <- TFs_data[TFs_data$Term==TF_David,'Genes']
        enriched_genes <- strsplit(enriched_genes, ', ')
        z<-data.frame('deseq' = enriched_genes[[1]])
        enriched_genes2<- as.character(z$deseq[z$deseq %in% DESeq$ENSEMBL])
        x <- as.numeric(exprs[exprs$symbol==TF_exprs,c(4:ncol(exprs))])
        cors <- c()
        
        random_genes <- sample(exprs$ENSEMBL, length(enriched_genes2)*4, replace=FALSE) # Here is random sampling
        for (gene in random_genes) {
          if (gene %in% exprs$ENSEMBL){
            y <- as.numeric(exprs[exprs$ENSEMBL==gene,c(4:ncol(exprs))])
            cors<- c(cors, cor.test(x,y, method='spearman')$estimate)}
        }
        
        return(median(na.omit(cors)))
      }
      
      result_random <- setNames(data.frame(matrix(ncol = 3, nrow = 0), stringsAsFactors=F), c("TF_David", "TF_exprs", "median_correlation"))
      
      for (tf in TFs_enriched) {
        if (tf %in% mapping$David_name){
          #print(tf)
          cors <- rysuj_random(tf)
          result_random <- rbind(result_random, data.frame('TF_David'=tf, 'TF_exprs'=mapping$exprs_name[mapping$David_name==tf], 'median_correlation'=cors))
        }
      }
      write.csv(arrange(result_random, desc(median_correlation)), paste0('./TADs/enriched_TFs/', reg_tad,'_TFs_correlations_random.csv'))
      
      
      ### ILLUSTRATE THE REAL AND RANDOM CORRELATIONS
      reg_tad <- '1101'
      random <- read.csv(paste0('./TADs/enriched_TFs/', reg_tad, '_TFs_correlations_random.csv'), header = T, stringsAsFactors = F, row.names = 1)
      real <- read.csv(paste0('./TADs/enriched_TFs/', reg_tad, '_TFs_correlations_onlyDE.csv'), header = T, stringsAsFactors = F, row.names = 1)
      data_to_plot <- real
      data_to_plot$random_median <-random$median_correlation[match(real$TF_David, random$TF_David)]
      height_to_move <- (max(random$median_correlation)-min(random$median_correlation))*0.2
      data_to_plot <- data.frame(data_to_plot[duplicated(data_to_plot$TF_exprs)==F,])

      ggplot(data_to_plot, aes(x=TF_exprs, y = median_correlation)) +
        geom_bar(stat = "identity") +
        scale_color_gradient2(low='cornflowerblue', high='coral', mid = 'grey', name='Spearman rho \n(TFs and \nrandom targets)')+
        geom_point(aes(x=TF_exprs, y=random_median, color=random_median, size=5))+
        xlab('Sorted TFs') +
        ylab('Median Spearman rho') +
        scale_x_discrete(limits=data_to_plot$TF_exprs) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
        geom_hline(yintercept = max(random$median_correlation), color = "coral", linetype = "dashed") +
        annotate("text", label = "Max random", x = length(data_to_plot$TF_exprs)/2, y = max(random$median_correlation)+height_to_move, color = "coral") +
        geom_hline(yintercept = min(random$median_correlation), color = "cornflowerblue", linetype = "dashed") +
        annotate("text", label = "Min random", x = length(data_to_plot$TF_exprs)/2, y = min(random$median_correlation)-height_to_move, color = "cornflowerblue")+
        guides(size=FALSE) +
        theme(text = element_text(size=15))
      ggsave(paste0('./TADs/enriched_TFs/', reg_tad,'_real_cors.png'), height =4, width=14)
      
      max_random <- max(random$median_correlation)
      min_random <- min(random$median_correlation)
      
      for (row in 1:nrow(result)) {
        tf <- result$TF_exprs[row]
        if (result$median_correlation[row] >0 ) {
          if (result$median_correlation[row] > max_random) {
            repeated_enriched_TFs <- rbind(repeated_enriched_TFs, data.frame('enriched_TF_above_random'=tf, 'TAD_nr'=reg_tad, 'median_correlation'=result$median_correlation[row]))
          }
        }
        
        if (result$median_correlation[row] <0 ) {
          if (result$median_correlation[row] < min_random) {
            repeated_enriched_TFs <- rbind(repeated_enriched_TFs, data.frame('enriched_TF_above_random'=tf, 'TAD_nr'=reg_tad, 'median_correlation'=result$median_correlation[row]))
          }
        }
      }
    }
    write.csv(repeated_enriched_TFs, './TADs/enriched_TFs/repeated_in_many_TADs_enriched_TFs_with_above_random_correlations.csv')
  }
}

repeated_enriched_TFs <- read.csv('./TADs/enriched_TFs/repeated_in_many_TADs_enriched_TFs_with_above_random_correlations.csv', stringsAsFactors = F, header=T)
table_rep_enr_TFs <- table(repeated_enriched_TFs$enriched_TF_above_random)
sort(table_rep_enr_TFs)
length(table_rep_enr_TFs)

# Corrplots: correlations of PCDHGA genes with TFs enriched for 1101 and confirmed by BMO:
BMO_genes <- c('BACH1', 'FOSL1', 'IRX3', 'ZBTB33', 'NFATC3', 'SOX1', 'SOX21', 'YY1', 'ZFP42', 'ZNF418', 'ZNF41') # KAISO = ZBTB33
PCDHG_genes <- c('ENSG00000081853', 'ENSG00000204956', 'ENSG00000254245', 'ENSG00000262576', 'ENSG00000253485', 'ENSG00000253953', 'ENSG00000253305', 
                 'ENSG00000276547', 'ENSG00000262209', 'ENSG00000253910', 'ENSG00000254221', 'ENSG00000253767', 'ENSG00000253537', 'ENSG00000253731')

PCDHG_genes_expressions <- data.frame(t(exprs[exprs$ENSEMBL %in% PCDHG_genes, c('symbol', colnames(exprs)[grepl('PA|DA|GB|PG', colnames(exprs))])]))
colnames(PCDHG_genes_expressions) <- PCDHG_genes_expressions[rownames(PCDHG_genes_expressions)=='symbol',]
PCDHG_genes_expressions <- PCDHG_genes_expressions[2:nrow(PCDHG_genes_expressions),] 
PCDHG_genes_expressions[,] %<>% lapply(function(x) as.numeric(as.character(x)))

BMO_TFs_genes_expressions <- data.frame(t(exprs[exprs$symbol %in% BMO_genes, c('symbol', colnames(exprs)[grepl('PA|DA|GB|PG', colnames(exprs))])]))
colnames(BMO_TFs_genes_expressions) <- BMO_TFs_genes_expressions[rownames(BMO_TFs_genes_expressions)=='symbol',]
BMO_TFs_genes_expressions <- BMO_TFs_genes_expressions[2:nrow(BMO_TFs_genes_expressions),] 
BMO_TFs_genes_expressions[,] %<>% lapply(function(x) as.numeric(as.character(x)))

corr_mat=cor(PCDHG_genes_expressions, BMO_TFs_genes_expressions, method="s")
cor_col = colorRampPalette(c("cornflowerblue", "red"))
png('corplot_for_PCDHG_genes_and_BMO_TFs.png')
corrplot(t(corr_mat), tl.col = 'black', col = cor_col(20))#, method='number')
dev.off()
