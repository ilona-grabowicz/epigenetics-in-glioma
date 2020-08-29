# Author: Ilona E. Grabowicz
# Date: 2019-2020

library(tidyr)
library(dplyr)
library(data.table)
library(ggrepel)
library(ggplot2)

######### Calculating correlations between gene expressions and epigenetic marks depositions on a TAD level
######### and visualizing the results.

dataset <- 'RNASeq'
comparison <-  'grade1_grade4'
part <- 'promoters'
##########################
funkcja <- function(dataset, comparison, part) {
  TADs <- read.csv('./input_files/TADborders-Fetal_brain_D_cortical_plates_Won_et_al_2016.bed', stringsAsFactors = F)
  genes <- read.csv(paste0('DeSeq2_', comparison, '.csv'), sep='\t')
  genes <- separate(genes,col = gene, into=c('ENSEMBL', 'symbol', 'location'), sep = '_') 
  genes <- separate(genes, col=location, into=c('chromosome', 'rest'), sep=':')
  genes <- separate(genes, col=rest, into=c('start', 'stop'))
  genes <- separate(genes, col=chromosome, into=c('chromosome', 'x'))
  genes$x<-NULL
  genes$start<-as.integer(genes$start)
  genes$stop<-as.integer(genes$stop)
  genes<- na.omit(genes[genes$FDR<0.01,])
  if (dataset == 'RNASeq'){
    DE <- genes }
  if (dataset!='RNASeq') {
    DE <- read.csv(paste0('DeSeq2_', dataset, '_', comparison,'_', part, '.csv'), sep='\t', stringsAsFactors = F)
    genes_coords <- genes[,c('ENSEMBL', 'symbol', 'chromosome', 'start', 'stop')]
    names(DE)[names(DE) == 'gene'] <- 'ENSEMBL'
    DE <- inner_join(x=DE, y=genes_coords, by = c("ENSEMBL" = "ENSEMBL"))
  }
  # Which genes are in which TADs?
  list_of_TADs <-list()
  TADs_with_symbols <- list()
  for (row in 1:nrow(TADs)) { #nrow(TADs)
    cat('.')
    my_DE <- DE[DE$chromosome == TADs[row,"chr"],]
    is.included <- mapply(function(de_start,de_stop,tad_start,tad_stop){de_start >= tad_start &  de_stop <= tad_stop}, my_DE$start, my_DE$stop, TADs[row,]$start, TADs[row,]$stop)
    list_of_TADs[[row]] <- my_DE[is.included,'log2_FC']
    TADs_with_symbols[[row]] <- my_DE[is.included,'ENSEMBL'] #symbol - then instead of ENSEMBL symbols are used
  }
  TADs_with_symbols <- lapply(TADs_with_symbols, function(x) {as.data.frame(t(x))})
  TADs_with_symbols <- lapply(TADs_with_symbols, function(x) { if(length(x)==0) data.frame(V1 = "xxx") else x} )
  TADs_with_symbols <- rbindlist(TADs_with_symbols, fill = T)
  TADs_with_symbols$V1[TADs_with_symbols$V1 == "xxx"] <- NA
  names(TADs_with_symbols) <- paste0("gene_",1:ncol(TADs_with_symbols))
  dim(TADs_with_symbols)
  means <- lapply(list_of_TADs, function(x) mean(x, na.rm = T)) #result is a list
  medians <- lapply(list_of_TADs, function(x) median(x, na.rm = T))
  TADs_df <- data.frame(TAD_nr=1:nrow(TADs), chr=TADs$chr, start=TADs$start, stop=TADs$stop, means=unlist(means), medians=unlist(medians), TADs_with_symbols)
  return(TADs_df)
}

rnaseq <- funkcja('RNASeq', 'grade1_grade23')
h3k4me3 <- funkcja('H3K4me3', 'grade1_grade23', 'promoters')
h3k27ac <- funkcja('H3K27ac', 'grade1_grade23', 'promoters')
regulated_TADs <- read.csv('./TADs/counting_DE_TADS.csv')
regulated_TADs <- regulated_TADs$Var1[regulated_TADs$Freq>=2]
rnaseq$regulated <- rnaseq$TAD_nr %in% regulated_TADs
# RNASeq vs H3K4me3:
data_for_plotting <- data.frame("TAD_nr" = rnaseq$TAD_nr, "RNA_means" = rnaseq$means, "h3k4me3" = h3k4me3$means, "regulated" = rnaseq$regulated)

# Checking data normality:
a<- rnorm(n=1000)
shapiro.test(a)
shapiro.test(data_for_plotting$RNA_means)
y <- rnorm(n=length(!is.na(data_for_plotting$RNA_means)), mean = mean(!is.na(data_for_plotting$RNA_means)), sd = sd(!is.na(data_for_plotting$RNA_means)))
qqplot(x=data_for_plotting$h3k4me3, y=y)
qqline(y)
# Distribution is not normal.

# Calculate correlation:
kor <- cor.test(data_for_plotting$RNA_means, data_for_plotting$h3k4me3, method='spearman')
kor
kor$p.value
kor <- cor.test(data_for_plotting$RNA_means[data_for_plotting$regulated==T], data_for_plotting$h3k4me3[data_for_plotting$regulated==T], method='spearman')
kor
kor$p.value
# Plot:
ggplot(data_for_plotting, aes(x=RNA_means, y=h3k4me3)) +
  geom_hline(yintercept=0, color = "grey") +
  geom_vline(xintercept=0, color = "grey") +
  geom_point(aes(color=regulated)) +
  scale_colour_manual(name="'glioma' \nTADs", values = c("TRUE" = "red","FALSE"="cornflowerblue")) +
  geom_point(data = data_for_plotting[data_for_plotting$regulated == TRUE,], col = "red") +
  geom_text_repel(data = data_for_plotting[data_for_plotting$regulated == TRUE,], label = data_for_plotting[data_for_plotting$regulated == TRUE, "TAD_nr"]) +
  xlab('Means of RNASeq log2(FC)') + ylab('Means of H3K4me3 \nlog2(FC) at promoters') +
  theme(text = element_text(size=15)) 
  #ggtitle('Means of TADs expressions and H3K4me3 \nhistone marks peaks,only DE genes')
ggsave('Means_RNAseq_h3k4me3_regulated_TADs_only_DE_genes_2orMore_grade1_23.png', width=5, height=3.5)

# RNASeq vs H3K27ac:
data_for_plotting <- data.frame("TAD_nr" = rnaseq$TAD_nr, "RNA_means" = rnaseq$means, "h3k27ac" = h3k27ac$means, "regulated" = rnaseq$regulated)

# Calculate correlation:
kor <- cor.test(data_for_plotting$RNA_means, data_for_plotting$h3k27ac, method='spearman')
kor
kor$p.value
kor <- cor.test(data_for_plotting$RNA_means[data_for_plotting$regulated==T], data_for_plotting$h3k27ac[data_for_plotting$regulated==T], method='spearman')
kor
kor$p.value

# Plot:
ggplot(data_for_plotting, aes(x=RNA_means, y=h3k27ac)) +
  geom_hline(yintercept=0, color = "grey") +
  geom_vline(xintercept=0, color = "grey") +
  geom_point(aes(color=regulated)) +
  scale_colour_manual(name="'glioma' \nTADs", values = c("TRUE" = "red","FALSE"="cornflowerblue")) +
  geom_point(data = data_for_plotting[data_for_plotting$regulated == TRUE,], col = "red") +
  geom_text_repel(data = data_for_plotting[data_for_plotting$regulated == TRUE,], label = data_for_plotting[data_for_plotting$regulated == TRUE, "TAD_nr"]) +
  xlab('Means of RNASeq log2(FC)') + ylab('Means of h3k27ac \nlog2(FC) at promoters') +
  theme(text = element_text(size=15)) 
  #ggtitle('Means of TADs expressions and h3k27ac \nhistone marks peaks, only DE genes')
ggsave('Means_RNAseq_h3k27ac_regulated_TADs_only_DE_genes_2orMore_grade1_23.png', width=5, height=3.5)

  