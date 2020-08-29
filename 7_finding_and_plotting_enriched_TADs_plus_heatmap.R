# Author: Ilona E. Grabowicz
# Date: 2019-2020 

library(tidyr)
library(dplyr)
library(data.table)
library(ggplot2)
library(scales)
library(ggrepel)

######## Finding the enriched glioma TADs and plotting a heatmap with them: 

#dataset <- 'H3K27me3'
#dataset <- 'RNASeq'
dataset <- 'ATACseq'
comparison <- 'grade1_grade23'
part <- 'promoters'
#############
# nr_of_genes_cutoff is a minimal number of genes within a TAD
funkcja <- function(dataset, comparison, part, nr_of_genes_cutoff) {
  nr_of_genes_cutoff<-3
  if (dataset=='RNASeq'){
    TADs_with_genes <- read.csv(paste0('./TADs/TADs_with_genes_ENSEMBL_and_FC_', dataset, '_', comparison, '.csv'), stringsAsFactors = F)
    DE_genes <- read.csv(paste0('./DeSeq2_', comparison, '.csv'), sep='\t')
    DE_genes <- separate(DE_genes,col = gene, into=c('ENSEMBL', 'symbol', 'location'), sep = '_') 
  }
  if (dataset!='RNASeq') {
    TADs_with_genes <- read.csv(paste0('./TADs/TADs_with_genes_ENSEMBL_and_FC_', dataset, '_', comparison, '_', part, '.csv'), stringsAsFactors = F)
    DE_genes <- read.csv(paste0('./DeSeq2_', dataset, '_', comparison, '_', part, '.csv'), sep='\t')
    names(DE_genes)[names(DE_genes) == 'gene'] <- 'ENSEMBL'
  }
  TADs_with_genes$X<- NULL
  DE_genes$significant <- DE_genes$FDR<0.01
  prob_in_genome <- length(which(DE_genes$significant==TRUE)) / nrow(DE_genes)
  percent <- lengths <- tad_number <- quantities <- probs <- c()
  #row<- 1
  for (row in 1:nrow(TADs_with_genes)) {#nrow(TADs_with_genes)
    cat('.')
    isna <- is.na(TADs_with_genes[row,grepl( "gene", names(TADs_with_genes)) ])
    gene_names <- TADs_with_genes[row,grepl( "gene", names(TADs_with_genes)) ][!isna]
    if (length(gene_names)>=nr_of_genes_cutoff) { # Dla TADow, co maja conajmniej x genow
      if_DE <- gene_names %in% DE_genes$ENSEMBL[DE_genes$significant==TRUE]
      tad_number <- c(tad_number, TADs_with_genes[row, "TAD_nr"])
      percent_single <- length(if_DE[if_DE==TRUE]) / (length(if_DE[if_DE==TRUE]) + length(if_DE[if_DE==FALSE]) ) *100 
      percent <- c(percent,percent_single)
      quantity <- length(if_DE[if_DE==TRUE])
      quantities <- c(quantities, quantity)
      dlugosc <-  length(gene_names)
      lengths <- c(lengths, dlugosc)
      probs <- c(probs, pbinom(quantity-1, dlugosc, p = prob_in_genome, lower.tail = F))
    }
  }
  new_table<- data.frame(tad_number, lengths, quantities, percent, probs)
  new_table$BH <-  p.adjust(new_table$probs, method = 'BH')
  new_table$significant <- 'insignificant'
  new_table$significant[new_table$BH<=0.05] <- 'significant'
  new_table <- arrange(new_table, probs)
  write.csv(new_table, paste0('./TADs/Significant_TADs_', dataset, '_', comparison,'_',part, '.csv'))
  # g <- ggplot(new_table, aes(lengths, percent, color=significant)) +
  #   geom_point() +
  #   ylab('Share of DE genes [%]') +
  #   xlab('Number of genes in a TAD') +
  #   scale_color_discrete(name="Enrichment") +
  #   #ggtitle('TADs enriched\nBinomial distribution test significance,\nBH correction p-val<=0.05') +
  #   geom_text_repel(data=new_table[new_table$significant=='significant'|new_table$significant=='very significant',], aes(label=tad_number), color=c('#00BFC4'))
  # print(g)
  # ggsave(paste0('./TADs/Number_of_genes_per_TAD_vs_percent_of_DE_genes_', dataset, '_', comparison,'_',part, '.png'), width = 4.5, height = 3)
  return(new_table$tad_number[new_table$BH<0.05])
}
rna_grade1_grade4 <- funkcja('RNASeq', 'grade1_grade4', 'promoters', 3)
rna_grade1_grade23 <- funkcja('RNASeq', 'grade1_grade23', 'promoters', 3)
#rna_normal_grade1 <- funkcja('RNASeq', 'normal_grade1', 'promoters', 3)
rna_grade23_grade4 <- funkcja('RNASeq', 'grade23_grade4', 'promoters', 3)

atac_grade1_grade4_prom <- funkcja('ATACseq', 'grade1_grade4', 'promoters', 3)
atac_grade1_grade23_prom <- funkcja('ATACseq', 'grade1_grade23', 'promoters', 3)
atac_grade23_grade4_prom <- funkcja('ATACseq', 'grade23_grade4', 'promoters', 3)

atac_grade1_grade4_genes <- funkcja('ATACseq', 'grade1_grade4', 'genes', 3)
atac_grade1_grade23_genes <- funkcja('ATACseq', 'grade1_grade23', 'genes', 3)
atac_grade23_grade4_genes <- funkcja('ATACseq', 'grade23_grade4', 'genes', 3)

dname_grade1_grade4_prom <- funkcja('DNAme', 'grade1_grade4', 'promoters', 3)
dname_grade1_grade23_prom <- funkcja('DNAme', 'grade1_grade23', 'promoters', 3)
dname_grade23_grade4_prom <- funkcja('DNAme', 'grade23_grade4', 'promoters', 3)

dname_grade1_grade4_genes <- funkcja('DNAme', 'grade1_grade4', 'genes', 3)
dname_grade1_grade23_genes <- funkcja('DNAme', 'grade1_grade23', 'genes', 3)
dname_grade23_grade4_genes <- funkcja('DNAme', 'grade23_grade4', 'genes', 3)

dnase_grade1_grade4_prom <- funkcja('DNase', 'grade1_grade4', 'promoters', 3)
dnase_grade1_grade23_prom <- funkcja('DNase', 'grade1_grade23', 'promoters', 3)
dnase_grade23_grade4_prom <- funkcja('DNase', 'grade23_grade4', 'promoters', 3)

dnase_grade1_grade4_genes <- funkcja('DNase', 'grade1_grade4', 'genes', 3)
dnase_grade1_grade23_genes <- funkcja('DNase', 'grade1_grade23', 'genes', 3)
dnase_grade23_grade4_genes <- funkcja('DNase', 'grade23_grade4', 'genes', 3)

h3k4me3_grade1_grade4_prom <- funkcja('H3K4me3', 'grade1_grade4', 'promoters', 3)
h3k4me3_grade1_grade23_prom <- funkcja('H3K4me3', 'grade1_grade23', 'promoters', 3)
h3k4me3_grade23_grade4_prom <- funkcja('H3K4me3', 'grade23_grade4', 'promoters', 3)

h3k27ac_grade1_grade4_prom <- funkcja('H3K27ac', 'grade1_grade4', 'promoters', 3)
h3k27ac_grade1_grade23_prom <- funkcja('H3K27ac', 'grade1_grade23', 'promoters', 3)
h3k27ac_grade23_grade4_prom <- funkcja('H3K27ac', 'grade23_grade4', 'promoters', 3)

h3k27ac_grade1_grade4_genes <- funkcja('H3K27ac', 'grade1_grade4', 'genes', 3)
h3k27ac_grade1_grade23_genes <- funkcja('H3K27ac', 'grade1_grade23', 'genes', 3)
h3k27ac_grade23_grade4_genes <- funkcja('H3K27ac', 'grade23_grade4', 'genes', 3)

all <- list('rna_grade1_grade4' = rna_grade1_grade4, 
            'rna_grade1_grade23' = rna_grade1_grade23, 
            'rna_grade23_grade4' = rna_grade23_grade4,
            'atac_grade1_grade4_prom' = atac_grade1_grade4_prom, 
            'atac_grade1_grade23_prom' = atac_grade1_grade23_prom, 
            'atac_grade23_grade4_prom' = atac_grade23_grade4_prom,
         'atac_grade1_grade4_genes' = atac_grade1_grade4_genes, 
         'atac_grade1_grade23_genes' = atac_grade1_grade23_genes, 
         'atac_grade23_grade4_genes' = atac_grade23_grade4_genes,
         'dname_grade1_grade4_prom' = dname_grade1_grade4_prom, 
         'dname_grade1_grade23_prom' = dname_grade1_grade23_prom, 
         'dname_grade23_grade4_prom' = dname_grade23_grade4_prom,
         'dname_grade1_grade4_genes' = dname_grade1_grade4_genes, 
         'dname_grade1_grade23_genes' = dname_grade1_grade23_genes, 
         'dname_grade23_grade4_genes' = dname_grade23_grade4_genes,
         'dnase_grade1_grade4_prom' = dnase_grade1_grade4_prom, 
         'dnase_grade1_grade23_prom' = dnase_grade1_grade23_prom, 
         'dnase_grade23_grade4_prom' = dnase_grade23_grade4_prom,
         'dnase_grade1_grade4_genes' = dnase_grade1_grade4_genes, 
         'dnase_grade1_grade23_genes' = dnase_grade1_grade23_genes, 
         'dnase_grade23_grade4_genes' = dnase_grade23_grade4_genes,
         'h3k4me3_grade1_grade4_prom' = h3k4me3_grade1_grade4_prom, 
         'h3k4me3_grade1_grade23_prom' = h3k4me3_grade1_grade23_prom, 
         'h3k4me3_grade23_grade4_prom' = h3k4me3_grade23_grade4_prom,
         'h3k27ac_grade1_grade4_prom' = h3k27ac_grade1_grade4_prom, 
         'h3k27ac_grade1_grade23_prom' = h3k27ac_grade1_grade23_prom, 
         'h3k27ac_grade23_grade4_prom' = h3k27ac_grade23_grade4_prom,
         'h3k27ac_grade1_grade4_genes' = h3k27ac_grade1_grade4_genes, 
         'h3k27ac_grade1_grade23_genes' = h3k27ac_grade1_grade23_genes, 
         'h3k27ac_grade23_grade4_genes' = h3k27ac_grade23_grade4_genes)
sort(table(unlist(all)))
write.csv(sort(table(unlist(all))), './TADs/counting_DE_TADS.csv')

length_function <- function(x) {
  vec <- c()
  for (el in x) {
    vec <- c(vec,length(el))
  }
  return(vec)
}
length_function(all)
max_length <- max(length_function(all))
df <- data.frame(matrix(nrow = max_length, ncol = length(all)))
for (el in 1:length(all)) {
  df[,el]<- c(all[[el]], rep(NA, max_length-length(all[[el]])))
  colnames(df)[el] <- names(all[el])
  }
write.csv(df , './TADs/TAD_enrichments.csv')

### Heatmap
length(colnames(df))/3 #[grepl('rna|h3k4me3|h3k27ac|dnase|dname|atac', colnames(df))]
experiments <- c('rna', 'h3k4me3', 'h3k27ac', 'dnase', 'dname', 'atac')
parts <- c('prom', 'genes')
grades <- c('grade1_grade23', 'grade23_grade4', 'grade1_grade4')
#grades <- c('grade1_grade23')
 grade <- 'grade1_grade4'
# part <- 'genes'
# exp <- 'rna'
 draw_heatmap <- function(){
   my_table<-data.frame(matrix(0, ncol= length(colnames(df))/3, nrow = length(unique(unlist(all)))))
   rownames(my_table) <- unique(unlist(all))
   colnames(my_table) <- c('rna', 'h3k4me3_prom', 'h3k27ac_genes', 'h3k27ac_prom', 'dnase_genes', 'dnase_prom', 
                           'dname_genes', 'dname_prom', 'atac_genes', 'atac_prom')
   for (grade in grades) {
     for (exp in experiments) {
       for (part in parts)  {
         tmp <- df [ , grepl(part, colnames(df)) & grepl(exp, colnames(df)) & grepl(grade, colnames(df))] #tads is df
         wzbogacone_tady <- tmp[!is.na(tmp)]
         my_table[rownames(my_table)%in%wzbogacone_tady , grepl(part, colnames(my_table)) & grepl(exp, colnames(my_table))] <- my_table[rownames(my_table)%in%wzbogacone_tady , grepl(part, colnames(my_table)) & grepl(exp, colnames(my_table))] +1
       }
       if (exp=='rna'){
         tmp <- df [ , grepl(exp, colnames(df)) & grepl(grade, colnames(df))] #tads is df
         wzbogacone_tady <- tmp[!is.na(tmp)]
         my_table[rownames(my_table)%in%wzbogacone_tady , grepl(exp, colnames(my_table))] <- my_table[rownames(my_table)%in%wzbogacone_tady , grepl(exp, colnames(my_table))] +1
       }
     }
   }
   my_table$tad <- rownames(my_table)
   ktore <- apply(my_table[,c(1:ncol(my_table)-1)], 1, sum)>=2
   my_table1 <- my_table[ktore,]
   colnames(my_table1) <- c('RNAseq', 'H3K4me3_prom', 'H3K27ac_genes', 'H3K27ac_prom', 'DNAseI_genes', 'DNAseI_prom', 
                            'DNAme_genes', 'DNAme_prom', 'ATACseq_genes', 'ATACseq_prom', 'tad')
   df_for_plot <- melt(my_table1,id.vars = c("tad"))
   sort(as.numeric(df_for_plot$tad))
   df_for_plot$tad <- factor(df_for_plot$tad,levels = c(sort(unique(as.numeric(df_for_plot$tad)))))
   class(df_for_plot$value)
   kolory <- c("#101C30","#6A7A9C","#AEC9FF", "#D9ECFF")
   p <- ggplot(df_for_plot,aes(x=tad,y=variable, fill=factor(value))) +
     geom_tile(colour="white",size=0.25) + labs(x = "TAD numbers", y = "Chromatin assay") +#, title = "TADs with exceptionally high proportion of \nDifferentially Active genes") +
     scale_y_discrete(expand=c(0,0)) +
     coord_fixed()+ #maintains aspect ratio = kwadraty
     scale_fill_manual(values = kolory ,na.value="grey90", name='Enriched between \nnumber of grades') +
     theme_grey(base_size=12)+ #text size
     theme(
       #bold font for both axis text
       axis.text=element_text(face="bold"),
       #set thickness of axis ticks
       axis.ticks=element_line(size=0.4),
       #remove plot background
       plot.background=element_blank(),
       #remove plot border
       panel.border=element_blank()) +
     theme(axis.text.x = element_text(angle = 45, hjust = 1))
   
   print(p)
   ggsave('./TADs/TADs_heatmap_all.png' , width=8, height = 2.9)
   return(levels(df_for_plot$tad))
 }
 kolejnosc_x_axis <- draw_heatmap()

###########################
# Plots for each grade separately:
###########################
for (grade in grades) {
  my_table<-data.frame(matrix(0, ncol= length(colnames(df))/3, nrow = length(kolejnosc_x_axis)))
  rownames(my_table) <- kolejnosc_x_axis
  colnames(my_table) <- c('rna', 'h3k4me3_prom', 'h3k27ac_genes', 'h3k27ac_prom', 'dnase_genes', 'dnase_prom', 
                          'dname_genes', 'dname_prom', 'atac_genes', 'atac_prom')
    for (exp in experiments) {
      for (part in parts)  {
        tmp <- df [ , grepl(part, colnames(df)) & grepl(exp, colnames(df)) & grepl(grade, colnames(df))] #tads is df
        wzbogacone_tady <- tmp[!is.na(tmp)]
        my_table[rownames(my_table)%in%wzbogacone_tady , grepl(part, colnames(my_table)) & grepl(exp, colnames(my_table))] <- my_table[rownames(my_table)%in%wzbogacone_tady , grepl(part, colnames(my_table)) & grepl(exp, colnames(my_table))] +1
      }
      if (exp=='rna'){
        tmp <- df [ , grepl(exp, colnames(df)) & grepl(grade, colnames(df))] #tads is df
        wzbogacone_tady <- tmp[!is.na(tmp)]
        my_table[rownames(my_table)%in%wzbogacone_tady , grepl(exp, colnames(my_table))] <- my_table[rownames(my_table)%in%wzbogacone_tady , grepl(exp, colnames(my_table))] +1
      }
    }
  
  my_table$tad <- rownames(my_table)
  my_table1 <- my_table
  colnames(my_table1) <- c('RNAseq', 'H3K4me3_prom', 'H3K27ac_genes', 'H3K27ac_prom', 'DNAseI_genes', 'DNAseI_prom', 
                           'DNAme_genes', 'DNAme_prom', 'ATACseq_genes', 'ATACseq_prom', 'tad')
  df_for_plot <- melt(my_table1,id.vars = c("tad"))
  sort(as.numeric(df_for_plot$tad))
  df_for_plot$tad <- factor(df_for_plot$tad,levels = kolejnosc_x_axis)
  class(df_for_plot$value)
  if (grade=='grade1_grade23') {kolory <- c("#101C30", "#97CC5E")
  grade_title <- 'PA vs DA'}
  if (grade=='grade23_grade4') {kolory <- c("#101C30", "#B99EFF") 
  grade_title <- 'DA vs GB/PG'}
  if (grade=='grade1_grade4') {kolory <- c("#101C30", "#FF794D") 
  grade_title <- 'PA vs GB/PG'}
  
  p <- ggplot(df_for_plot,aes(x=tad,y=variable, fill=factor(value))) +
    geom_tile(colour="white",size=0.25) + labs(x = "TAD numbers", y = "Chromatin assay", title = grade_title) +
    scale_y_discrete(expand=c(0,0)) +
    coord_fixed()+ #maintains aspect ratio = squares
    scale_fill_manual(values = kolory , labels = NULL, na.value="grey90") +
    guides(fill = 'none') +
    theme_grey(base_size=15)+ 
    theme(
      #bold font for both axis text
      axis.text=element_text(face="bold"),
      #set thickness of axis ticks
      axis.ticks=element_line(size=0.4),
      #remove plot background
      plot.background=element_blank(),
      #remove plot border
      panel.border=element_blank()) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p)
  ggsave(paste0('./TADs/TADs_heatmap_', grade, '.png' ), width=7, height = 3)
}


