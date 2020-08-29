# Author: Ilona E. Grabowicz
# Date: 2019-2020
library(tidyr)
library(dplyr)
library(data.table)
library(ggplot2)
library(ggrepel)

# Script with analyses and visualisations of contacts of genes with enhancers, enhancers acetylations and DNA methylations.

### Matching DE genes with predicted long-range contacts with enhancers:
de_1_4_orig <- read.csv('DeSeq2_grade1_grade4.csv', stringsAsFactors = F, header = T, sep='\t')
de_1_4_orig <- separate(de_1_4_orig,col = gene, into=c('gene', 'symbol', 'loc'), sep = '_')
de_1_4_orig <- separate(de_1_4_orig,col = loc, into=c('chr', 'after'), sep = ':')
de_1_4_orig <- separate(de_1_4_orig,col = after, into=c('start', 'stop'), sep = '-')

de_1_4_orig$chr <- substr(de_1_4_orig$chr, 1, nchar(de_1_4_orig$chr)-1)
de_1_4 <- na.omit(de_1_4_orig[de_1_4_orig$FDR<=0.01,])

files = list.files(pattern="*.gff", path = './long_range_interactions/')
kontakty = data.frame(do.call(rbind, lapply(paste0('./long_range_interactions/', files), fread, skip=1)))
kontakty <- kontakty[,c(1,3,4,5,9)]
colnames(kontakty) <- c('chr', 'type_of_region', 'start', 'stop', 'bin')
kontakty <- separate(kontakty, col=bin, into=c('remove', 'bin'), sep='=')
kontakty$remove<-NULL
kontakty_predictions <- kontakty[kontakty$type_of_region=='prediction',] # Narrowing down to regions which would be in genes. ('region' means enhancer)
enhancery <- read.csv('./long_range_interactions/H3K27ac_coverage_normalised_on_H3K27ac_peaks_outside_promoters.tsv', header=T, sep='\t', stringsAsFactors = F)
enhancery$X <- NULL
wynik <- setNames(data.frame(matrix(ncol = 10, nrow = 0)), 
                  c('DE_ENSEMBL', 'symbol', 'chr', 'start', 'stop', 'prediction_start', 'prediction_stop', 'bin_ID', 'region_start', 'region_stop'))

i<-1
j<-1
for (i in 1:nrow(de_1_4)) { # 1:nrow(de_1_4)
  de_gene <- de_1_4$gene[i]
  de_symbol <- de_1_4$symbol[i]
  chrom <- de_1_4$chr[i]
  if (chrom=='chrY'){
    next()
  }
  start <- as.numeric(de_1_4$start[i])
  stop <- as.numeric(de_1_4$stop[i])
  tmp.kontakty <-kontakty[kontakty$chr== chrom,] # zawezam do tego samego chromosomu co jest gen w petli
  tmp.kontakty_pred <-tmp.kontakty[tmp.kontakty$type_of_region== 'prediction',] # zawezam tylko do przewidzianych miejsc (predykcji)
  for (j in 1:nrow(tmp.kontakty_pred)) { # nrow(tmp.kontakty_pred) # Iteruje przez predykcje
    #print(paste0('i',i, 'j',j))
    kont_chrom <- tmp.kontakty_pred$chr[j] # to w zasadzie nie powinno byc potrzebne?
    kont_start <- as.numeric(tmp.kontakty_pred$start[j])
    kont_stop <- as.numeric(tmp.kontakty_pred$stop[j])

    if (chrom == kont_chrom){
        if (  !(kont_stop<start | kont_start>stop)  ){ # tutaj znajduje predykcje ktora ma overlap z genem
          bin <- tmp.kontakty_pred$bin[j]
          region_start <- kontakty$start[kontakty$type_of_region=='region'&kontakty$bin==bin] # tj koordynata regionu (enhancera), ktory sie kontaktuje z predycja
          region_stop <- kontakty$stop[kontakty$type_of_region=='region'&kontakty$bin==bin]
          #na raz zrobienie df
          for (region_nr in 1:length(region_start)){ # czasami jeste wiecej regionow co ma polaczenie z ta sama predykcja, stad rozbijam je na pojedyncze wiersze
            tmp <- c(de_gene, de_symbol, chrom, start, stop, kont_start, kont_stop, bin, region_start[region_nr], region_stop[region_nr])
            wynik[nrow(wynik)+1,] <- tmp #Adding found matches to the result list 
          }
        }
      }
    }
  }

write.csv(wynik, 'DE_genes_and_contact_predictions.csv')
wynik$region_start<-as.integer(wynik$region_start)
wynik$region_stop<-as.integer(wynik$region_stop)

# Counting how many enhancer have predicted contacts with DE genes:
kont_reg<-kontakty[kontakty$type_of_region=='region',]
write.csv(kont_reg, 'Regiony_co_sa_w_kontakcie_z_predykacjami_overlapujacymi_z_DEgenes.csv') 

length(unique(kont_reg$start)) # That many enhancers ('regions') have predicted contact 
length(unique(kont_reg$start)) / nrow(enhancery) *100 # That percentage of enhancers have predicted contact
length(unique(kont_reg$stop)) # almost the same as kont_reg$start

length(unique(wynik$region_start)) # That many enhancers ('regions') have predicted contact with DE genes
length(unique(wynik$region_stop))
length(unique(wynik$region_start)) / nrow(enhancery) *100 # # That percentage of enhancers have predicted contact with DE genes
reg_DE <- unique(wynik$region_start) # Number of unique enhancers - For drawing a Venn diagram
  
# Counting how many DE genes have contacts with enhancers 

length(unique(wynik$DE_ENSEMBL)) # So many unique DE genes have predicted contacts with enhancers = 1716
length(unique(wynik$DE_ENSEMBL)) / nrow(de_1_4) *100 # 41% Such percentage of DE genes -||-


# Comparing expression levels of genes contacting enhancers vs not-contacting enhancers
wilcox.test(de_1_4$mean_both[de_1_4$gene%in%wynik$DE_ENSEMBL], # in-contact
              de_1_4$mean_both[!de_1_4$gene%in%wynik$DE_ENSEMBL])$p.value #not in contact
# How much?

a <- mean(de_1_4$mean_both[de_1_4$gene%in%wynik$DE_ENSEMBL]) 
b <- mean(de_1_4$mean_both[!de_1_4$gene%in%wynik$DE_ENSEMBL])
(a-b)/b

# Finding how many enhancers are differentially acetylated (H3K27ac)

wilcox <- c()
for (e in 1:nrow(enhancery)) { #nrow(enhancery) # file with enhancer counts
a <- (wilcox.test(as.numeric(enhancery[e,grepl('PA', colnames(enhancery))]), as.numeric(enhancery[e,grepl('GB|PG', colnames(enhancery))])))
wilcox <- c(wilcox, a$p.value)
}
enhancery$wilcox_pvalue <- wilcox

# Finding differentially acetylated enhancers being in contact with DE genes?
wynik$enhancer_pvalue <- wynik$enhancer_stop <- wynik$enhancer_start <- NA
l<-1
for (l in 1:nrow(wynik)) { #1:nrow(wynik) iterating through regions (enhancers) which have contact with 'predictions' that have overlap with genes
  chromosom <- wynik$chr[l]
  enhancery_sub <- enhancery[enhancery$chr==chromosom,]
  if (wynik$region_start[l]%in%enhancery_sub$start){
  wynik$enhancer_start[l] <- enhancery_sub$start[enhancery_sub$start==wynik$region_start[l]]
  wynik$enhancer_stop[l] <- enhancery_sub$end[enhancery_sub$start==wynik$region_start[l]]
  wynik$enhancer_pvalue[l] <- enhancery_sub$wilcox_pvalue[enhancery_sub$start==wynik$region_start[l]]
  }
}


wynik <- arrange(wynik, enhancer_pvalue)
write.csv(wynik, 'DE_genes_contact_predictions_and_enhancers_acetylations.csv')
#wynik <- read.csv('DE_genes_contact_predictions_and_enhancers_acetylations.csv')
#gene <- wynik$DE_ENSEMBL[wynik$enhancer_pvalue<=0.01]

# How many enhancers are differentially acetylated?
diff_acet <- nrow(enhancery[enhancery$wilcox_pvalue<0.01,]) 
diff_acet # That many enhancers were found to have a differential acetylation of H3K27

# That many enhancers were found to have a differential acetylation of H3K27 and are in contact with DE genes:
diff_a_reg <- length(unique(wynik$bin_ID[!is.na(wynik$enhancer_pvalue)&wynik$enhancer_pvalue<0.01&!is.na(wynik$region_start)] )) 
diff_a_reg / diff_acet *100 # Share of diff. acet. enhancers contacting DE genes per all diff. acet enhancers
diff_acet_starts <- enhancery$start[enhancery$wilcox_pvalue<0.01] #Start locations of diff. acet. enhancers

# Counting how many DE genes are in contact with differentially acetylated enhancers
length(unique(wynik$DE_ENSEMBL)) 
DE_enh <- c(unique(wynik$DE_ENSEMBL))
length(unique(wynik$DE_ENSEMBL)) / nrow(de_1_4) *100 # 41%
length(unique(wynik$DE_ENSEMBL[!is.na(wynik$enhancer_pvalue) & wynik$enhancer_pvalue<0.01])) # 117 That many DE genes have contact with differentially acetylated enhancer
DE_diff_acet_enh <- c(unique(wynik$DE_ENSEMBL[!is.na(wynik$enhancer_pvalue) & wynik$enhancer_pvalue<0.01])) # DE genes in contact wtit diff. acet. enhancers
  
# Venn enhancers
write.csv(reg_DE, 'reg_DE.csv',row.names=FALSE) # Number of unique enhancers
write.csv(diff_acet_starts, 'diff_acet_starts.csv',row.names=FALSE)
system("python overlap_between_enh_in_contact_with_DE_and_DiffAcetEnh.py")

#### Illustrating DE genes that have multiple contacts
# DE genes with contacts:
tail(sort(table(wynik$symbol[!is.na(wynik$enhancer_start)])))
geny_kontakty <- data.frame(table(wynik$symbol[!is.na(wynik$enhancer_start)])) 
geny_kontakty <- arrange(geny_kontakty, Freq)
geny_kontakty$Var1 <- as.character(geny_kontakty$Var1)
geny_kontakty$l_porz <- 1:nrow(geny_kontakty)
write.csv(geny_kontakty, 'geny_kontakty_grade1_4.txt', row.names=FALSE)

# Correlate number of contacts with mean exprs level: 
geny_kontakty_1_4 <- read.csv('geny_kontakty_grade1_4.txt', header=T, stringsAsFactors = F)
#geny_kontakty_1_23 <- read.csv('geny_kontakty_grade1_23.txt', header=T, stringsAsFactors = F)
#geny_kontakty_23_4 <- read.csv('geny_kontakty_grade23_4.txt', header=T, stringsAsFactors = F)
geny_kontakty_corr <- geny_kontakty_23_4
rna_original_mean <- rna_original

colnames(geny_kontakty_corr) <- c('symbol', 'Freq', 'l_porz')
rna_original_mean <- rna_original_mean[rna_original_mean$symbol%in%DE_genes_in_any_comparison,]
rna_original_mean$mean_PA <- rowMeans(subset(rna_original_mean, select = c(grepl('PA', colnames(rna_original_mean)))), na.rm = TRUE )
rna_original_mean$mean_DA <- rowMeans(subset(rna_original_mean, select = c(grepl('DA', colnames(rna_original_mean)))), na.rm = TRUE )
rna_original_mean$mean_GB <- rowMeans(subset(rna_original_mean, select = c(grepl('GB|PG', colnames(rna_original_mean)))), na.rm = TRUE )

geny_kontakty_corr$mean_expr_PA<-rna_original_mean$mean_PA[match(geny_kontakty_corr$symbol, rna_original_mean$symbol)]
geny_kontakty_corr$mean_expr_DA<-rna_original_mean$mean_DA[match(geny_kontakty_corr$symbol, rna_original_mean$symbol)]
geny_kontakty_corr$mean_expr_GB<-rna_original_mean$mean_GB[match(geny_kontakty_corr$symbol, rna_original_mean$symbol)]

geny_kontakty_corr$Freq_bin[geny_kontakty_corr$Freq==1] <- '1'
geny_kontakty_corr$Freq_bin[geny_kontakty_corr$Freq>1 & geny_kontakty_corr$Freq<=5] <- '1-5'
geny_kontakty_corr$Freq_bin[geny_kontakty_corr$Freq>5 & geny_kontakty_corr$Freq<=10] <- '6-10'
geny_kontakty_corr$Freq_bin[geny_kontakty_corr$Freq>10 & geny_kontakty_corr$Freq<=20] <- '11-20'
geny_kontakty_corr$Freq_bin[geny_kontakty_corr$Freq>20] <- '>20'
geny_kontakty_corr$Freq <- geny_kontakty_corr$l_porz <- NULL
geny_kontakty_corr <- melt(geny_kontakty_corr, by=c('symbol', 'Freq_bin'))
geny_kontakty_corr$Freq_bin <- factor(geny_kontakty_corr$Freq_bin, levels = c("1", "1-5", "6-10", "11-20", ">20"))
ggplot(geny_kontakty_corr, aes(x=Freq_bin, y=value, fill=variable)) +
  geom_boxplot() +
  scale_y_sqrt(limits=c(0, 10000))+
  xlab('Number of contacts') +
  ylab('Mean expression level') +
  scale_fill_manual(values = c('#56B4E9','#E69F00','#999999'), name=' ', labels = c('PA', 'DA', 'GB/PG')) +
  theme(text = element_text(size=15))

ggsave('Boxplots_comparing_contactNumbers_with_DE_Expression_23_4.png', width=5, height=4.5)
w <- wilcox.test(geny_kontakty_corr$value[geny_kontakty_corr$Freq_bin=='1' | geny_kontakty_corr$Freq_bin=='1-5'], 
            geny_kontakty_corr$value[geny_kontakty_corr$Freq_bin=='6-10' | geny_kontakty_corr$Freq_bin=='11-20' | geny_kontakty_corr$Freq_bin=='>20'])
w$p.value

# Number of enhancers contacting PRDM16 gene:
sum(!is.na(unique(wynik$enhancer_start[wynik$symbol=='PRDM16'])))

# DE genes with at least 2 contacts:
geny_kontakty_2 <- geny_kontakty[geny_kontakty$Freq>=2,]
write.csv(geny_kontakty_2, 'geny_kontakty2.csv')
ggplot(geny_kontakty_2, aes(x=l_porz, y = Freq)) +
  geom_point() +
  geom_label_repel(data = geny_kontakty_2[geny_kontakty_2$Freq>20,], #4 
                   label = geny_kontakty_2$Var1[geny_kontakty_2$Freq>20],
                   min.segment.length = 0, fontface = "italic", box.padding = 0.3) +
  xlab('Multi-loop genes') +
  ylab('Number of contacts') +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size=15))
ggsave('DE_genes_with_at_least2_contacts_grade1_4.png', width=4.5, height=4.5)
tail(geny_kontakty_2)
# Calculate how many DE genes have predicted contacts with multiple enhancers / multiple contacts with enhancer(s)
nrow(geny_kontakty_2) / nrow(de_1_4) *100


# Illustrating correlation between expression and acetylation of enhancer/enhancers or DNA methylation:

data <- read.csv('./input_files/RNAseq/combined_RNAs_normalised.csv', sep='\t', stringsAsFactors = F)#, col.names = 1)
data<- separate(data,col = gene, into=c('ENSEMBL', 'symbol', 'location'), sep = '_')
genes_in_contact_with_diffMeth_diffAcet_enhancers <- c('EGFR')#'GLI3', 'ITGB3BP', 'DNAH11', 'ADGRA1', 'ADGRF5P1', 'EGFR', 'EDN', 'FARSB', 'PHACTR1')
data_to_plot <- data.frame(matrix(nrow=0, ncol=3))
for (gen_tmp in genes_in_contact_with_diffMeth_diffAcet_enhancers) {
  pacjenci <- colnames(enhancery)[grepl('PA|DA|GB|PG', colnames(enhancery))]
  gen_tmp_exprs <- data[data$symbol==gen_tmp, pacjenci] # ekspresja tego genu
  
  enhancery_z_genem <- wynik[wynik$symbol==gen_tmp,]
  enh_tmp <- enhancery_z_genem[1,]
  enh_tmp_acet <- enhancery[enhancery$chr == enh_tmp$chr & enhancery$start == enh_tmp$region_start, pacjenci]
  data_to_plot_tmp <- inner_join(melt(enh_tmp_acet, variable.name = "patients", value.name = "enhancer acetylation"), 
                                 melt(gen_tmp_exprs, variable.name = "patients", value.name = "gene expression"), by='patients')
  data_to_plot <- rbind(data_to_plot, data_to_plot_tmp)
}
data_to_plot$grade[grepl('PA', data_to_plot$patients)] <- 'PA'
data_to_plot$grade[grepl('DA', data_to_plot$patients)] <- 'DA'
data_to_plot$grade[grepl('PG|GB', data_to_plot$patients)] <- 'GB/PG'
data_to_plot$grade <- factor(data_to_plot$grade, levels=c("PA", "DA", "GB/PG"))
cor <- cor.test(data_to_plot$`enhancer acetylation`, data_to_plot$`gene expression`, method='spearman')
ggplot(data_to_plot, aes(x = `enhancer acetylation` , y = `gene expression`, color = grade)) +
  geom_point(size=5, alpha=0.8) +
  theme(text = element_text(size=15),
        legend.title = element_blank()) +
  ggtitle(gen_tmp) +
  scale_color_manual(values=c('#56B4E9','#E69F00','#999999')) +
  scale_y_sqrt(labels=function(x) format(x, big.mark = " ", scientific = FALSE)) +
  scale_x_sqrt() +
  annotate("text", x = 10, y = 50000, label = paste0('rho = ', round(cor$estimate,2), # for EGFR: x = 10, y = 50000; for DNAH x = 10, y = 1000
                                                 '\n p = ', signif(cor$p.value, digits=3)))
ggsave(paste0(gen_tmp, '_geneExpression_enhancerAcetylation.png'))

# As above - with DNA methylation:
data <- read.csv('./input_files/RNAseq/combined_RNAs_normalised.csv', sep='\t', stringsAsFactors = F)#, col.names = 1)
data<- separate(data,col = gene, into=c('ENSEMBL', 'symbol', 'location'), sep = '_')
enhancer_methylations <- read.csv('./long_range_interactions/enhancer_methylations/enhancers_methylations.txt', header=T, sep='\t', stringsAsFactors = F)
genes_in_contact_with_diffMeth_diffAcet_enhancers <- c('EGFR')#'GLI3', 'ITGB3BP', 'DNAH11', 'ADGRA1', 'ADGRF5P1', 'EGFR', 'EDN', 'FARSB', 'PHACTR1')
data_to_plot <- data.frame(matrix(nrow=0, ncol=3))
for (gen_tmp in genes_in_contact_with_diffMeth_diffAcet_enhancers) {
  pacjenci <- colnames(enhancer_methylations)[grepl('PA|DA|GB|PG', colnames(enhancer_methylations))]
  gen_tmp_exprs <- data[data$symbol==gen_tmp, pacjenci] # ekspresja tego genu
  enhancery_z_genem <- wynik[wynik$symbol==gen_tmp,]
  enh_tmp <- enhancery_z_genem[1,]
  enh_tmp_meth <- enhancer_methylations[enhancer_methylations$chr == enh_tmp$chr & enhancer_methylations$enh_start == enh_tmp$region_start, pacjenci]
  data_to_plot_tmp <- inner_join(melt(enh_tmp_meth, variable.name = "patients", value.name = "Enhancer DNA methylation"), 
                                 melt(gen_tmp_exprs, variable.name = "patients", value.name = "Gene expression"), by='patients')
  data_to_plot <- rbind(data_to_plot, data_to_plot_tmp)
}
data_to_plot$grade[grepl('PA', data_to_plot$patients)] <- 'PA'
data_to_plot$grade[grepl('DA', data_to_plot$patients)] <- 'DA'
data_to_plot$grade[grepl('PG|GB', data_to_plot$patients)] <- 'GB/PG'
data_to_plot$grade <- factor(data_to_plot$grade, levels=c("PA", "DA", "GB/PG"))
data_to_plot <- na.omit(data_to_plot)
cor <- cor.test(data_to_plot$`Enhancer DNA methylation`, data_to_plot$`Gene expression`, method='spearman')
ggplot(data_to_plot, aes(x = `Enhancer DNA methylation` , y = `Gene expression`, color = grade)) +
  geom_point(size=5, alpha=0.8) +
  theme(text = element_text(size=15),
        legend.title = element_blank()) +
  ggtitle(gen_tmp) +
  scale_color_manual(values=c('#56B4E9','#E69F00','#999999')) +
  annotate("text", x = 0.6, y = 9000, label = paste0('rho = ', round(cor$estimate,2), # for EGFR x = 0.6, y = 9000; for DNAH x = 0.4, y = 1000
                                                     '\n p = ', signif(cor$p.value, digits=1)))
ggsave(paste0(gen_tmp, '_geneExpression_enhancerDNAMethylation.png'))

#### Illustrating DE genes that have contacts with differentially acetylated enhancers:
# DE genes in contact with diff. acet. enhancers:
geny_kontakty_rozn <- data.frame(table(wynik$symbol[!is.na(wynik$enhancer_pvalue) & wynik$enhancer_pvalue<0.01]))
geny_kontakty_rozn <- arrange(geny_kontakty_rozn, Freq)
geny_kontakty_rozn$Var1 <- as.character(geny_kontakty_rozn$Var1)
geny_kontakty_rozn$l_porz <- 1:nrow(geny_kontakty_rozn)
write.csv(geny_kontakty_rozn, 'geny_kontakty_rozn.txt')

ggplot(geny_kontakty_rozn, aes(x=l_porz, y = Freq)) +
  geom_point() +
  geom_label_repel(data = geny_kontakty_rozn[geny_kontakty_rozn$Freq>1,], 
                   label = geny_kontakty_rozn$Var1[geny_kontakty_rozn$Freq>1],
                   min.segment.length = 0,
                   max.iter = 10000,
                   fontface = "italic") +
  xlab('DE genes in contact(s) with \ndifferentially acetylated enhancers') +
  ylab('Number of contacts') +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size=15))
ggsave('DE_genes_with_diff_acetylated_enhancers.png', width=4.5, height=4.5)

# With at least 2 contacts:
geny_kontakty_rozn_2 <- geny_kontakty_rozn[geny_kontakty_rozn$Freq>=2,]
write.csv(geny_kontakty_rozn_2, 'geny_kontakty_rozn_2.txt')

## Illustrating enhancer acetylation with regard to WHO subtype:

PCDHGA_genes <- c('ENSG00000081853', 'ENSG00000204956', 'ENSG00000254245', 'ENSG00000262576', 'ENSG00000253485', 'ENSG00000253953', 'ENSG00000253305', 
                  'ENSG00000276547', 'ENSG00000262209', 'ENSG00000253910', 'ENSG00000254221', 'ENSG00000253767', 'ENSG00000253537', 'ENSG00000253731')

inv_gene <- 'ENSG00000081853'
enh_start <- wynik$region_start[wynik$DE_ENSEMBL==inv_gene & wynik$enhancer_pvalue<0.01 & !is.na(wynik$enhancer_pvalue)]
enh_start <- enh_start[1]
enh_stop <- wynik$region_stop[wynik$DE_ENSEMBL==inv_gene & wynik$enhancer_pvalue<0.01 & !is.na(wynik$enhancer_pvalue)]
enh_stop <- enh_stop[1]
enhancery[enhancery$start==enh_start & enhancery$end==enh_stop+1,]

data_to_plot <- enhancery[enhancery$start==enh_start & enhancery$end==enh_stop+1,]
data_to_plot1 <- melt(data_to_plot[,grepl('PA|DA|GB|PG', colnames(data_to_plot))])
data_to_plot1$tumor[grepl('PA',data_to_plot1$variable)] <- 'PA' 
data_to_plot1$tumor[grepl('DA',data_to_plot1$variable)] <- 'DA' 
data_to_plot1$tumor[grepl('GB',data_to_plot1$variable)] <- 'GB' 
data_to_plot1$tumor[grepl('PG',data_to_plot1$variable)] <- 'PG' 
data_to_plot1$tumor <- factor(data_to_plot1$tumor, levels=c("PA", "DA", "GB", "PG"))
ggplot(data_to_plot1, aes(x=tumor, y=value, colour = tumor)) +
  geom_point(size = 7, alpha=0.6) +
  xlab('Glioma type') + ylab('Enhancer acetylation coverage') +
  scale_x_discrete(limits=c('PA', 'DA', 'GB', 'PG')) +
  scale_color_manual(values = c('#56B4E9','#E69F00','#999999', 'black')) +
  theme(text = element_text(size=15),
        legend.title = element_blank())
  #ggtitle(paste0(wynik$symbol[wynik$DE_ENSEMBL==inv_gene],'\n', inv_gene, ' ','\nenhancer: ',data_to_plot$chr, ' ', data_to_plot$start, '-', data_to_plot$end))
ggsave(paste0(inv_gene,'_enhancer_acetylation.png'))
print(wilcox.test(as.numeric(data_to_plot1$value[grepl('PA',data_to_plot1$tumor)]) , as.numeric(data_to_plot1$value[grepl('GB|PG',data_to_plot1$tumor)])))

# Illustrating gene expression with regard to WHO subtype:

data <- read.csv('./input_files/RNAseq/combined_RNAs_normalised.csv', sep='\t', stringsAsFactors = F)#, col.names = 1)
data<- separate(data,col = gene, into=c('ENSEMBL', 'symbol', 'location'), sep = '_')
selected_genes <- inv_gene
for (which_gene in selected_genes){
  ktore <- which(data$ENSEMBL==which_gene)
  data1 <- data[ktore,]
  data1$location <- NULL
  expression <- melt(data1, id=c("ENSEMBL", "symbol"), variable.name = 'sample', value.name = 'expression', variable.factor=T)
  expression[grepl('PA',expression$sample)==TRUE, 'tumor'] <-'PA'
  expression[grepl('DA',expression$sample)==TRUE, 'tumor'] <-'DA'
  expression[grepl('GB',expression$sample)==TRUE, 'tumor'] <-'GB'
  expression[grepl('PG',expression$sample)==TRUE, 'tumor'] <-'PG'
  expression$tumor <- factor(expression$tumor, levels=c("PA", "DA", "GB", "PG"))
  ensembl <- unique(as.character(expression$ENSEMBL))
  symbol <- unique(as.character(expression$symbol))
  
  ggplot(expression, aes(x=tumor, y=expression, colour = tumor)) +
    geom_point(size = 7, alpha=0.6) +
    xlab('Glioma type') + ylab('Expression counts normalized') +
    scale_x_discrete(limits=c('PA', 'DA', 'GB', 'PG')) +
    scale_color_manual(values = c('#56B4E9','#E69F00','#999999', 'black')) +
    theme(text = element_text(size=15),
          legend.title = element_blank())
    ggsave(paste0(ensembl,'.png'))
  }

# Corrplot for all PCDHG genes expression and contacting enhancers acetylation:
library('corrplot')
library(magrittr)
#pacjenci <- colnames(enhancer_methylations)[grepl('PA|DA|GB|PG', colnames(enhancer_methylations))]
PCDHG_genes <- unique(wynik$symbol[!is.na(wynik$enhancer_start) & grepl('PCDHG', wynik$symbol)])
pacjenci <- pacjenci[pacjenci %in% colnames(enhancery)]
PCDHG_genes_expressions <- data.frame(t(data[data$symbol %in% PCDHG_genes, c('symbol', pacjenci)]))
colnames(PCDHG_genes_expressions) <- PCDHG_genes_expressions[rownames(PCDHG_genes_expressions)=='symbol',]
PCDHG_genes_expressions <- PCDHG_genes_expressions[2:nrow(PCDHG_genes_expressions),] 
PCDHG_genes_expressions[,] %<>% lapply(function(x) as.numeric(as.character(x)))
PCDHG_genes_enhancers_region_start <- unique(wynik$region_start[!is.na(wynik$enhancer_start) & grepl('PCDHG', wynik$symbol)])
PCDHG_genes_enhancers_region_chr <- unique(wynik$chr[!is.na(wynik$enhancer_start) & grepl('PCDHG', wynik$symbol)])
PCDHG_genes_enhancers <- data.frame(t(enhancery[enhancery$chr==PCDHG_genes_enhancers_region_chr & 
                                     enhancery$start %in% PCDHG_genes_enhancers_region_start, c( 'start', 'end', pacjenci)]))
colnames(PCDHG_genes_enhancers) <- paste0(PCDHG_genes_enhancers_region_chr, '_', 
                                          PCDHG_genes_enhancers[rownames(PCDHG_genes_enhancers)=='start',], '_', 
                                          PCDHG_genes_enhancers[rownames(PCDHG_genes_enhancers)=='end',])
PCDHG_genes_enhancers <- PCDHG_genes_enhancers[3:nrow(PCDHG_genes_enhancers),] 
PCDHG_genes_enhancers[,] %<>% lapply(function(x) as.numeric(as.character(x)))
cor_for_plot <- t(cor(PCDHG_genes_expressions, PCDHG_genes_enhancers))
cor_col = colorRampPalette(c("cornflowerblue", "red"))
png('corplot_for_PCHG_genes_and_3_enhancers.png')
corrplot(cor_for_plot, tl.col = 'black', col = cor_col(20))
dev.off()

# Correlating the expression level with enhancer acetylation of one of the PCHGA genes:
patients <- colnames(data1)[grepl('PA|DA|GB|PG',colnames(data1))]
patients <- patients[patients %in% colnames(data_to_plot) ]

cor.test(as.numeric(data1[,patients]), as.numeric(data_to_plot[,patients]), method='spearman')

# Correlation of the expression levels with enhancer acetylations for all the genes that have contacts with differentially acetylated enhancers:
data <- read.csv('./input_files/RNAseq/combined_RNAs_normalised.csv', sep='\t', stringsAsFactors = F)#, col.names = 1)
data <- separate(data,col = gene, into=c('ENSEMBL', 'symbol', 'location'), sep = '_')
data <- data[data$symbol%in%geny_kontakty_rozn$Var1,]

patients <- colnames(data)[grepl('PA|DA|GB|PG',colnames(data))]
patients <- patients[patients %in% colnames(enhancery) ]

kor_R <- kor_p <-c()
x <- y <- c()
for (gene in unique(geny_kontakty_rozn$Var1)) {
  inv_gene <- gene
  enh_start <- wynik$enhancer_start[wynik$symbol==inv_gene & wynik$enhancer_pvalue<0.01 & !is.na(wynik$enhancer_pvalue)]
  enh_start <- enh_start[1]
  enh_stop <- wynik$enhancer_stop[wynik$symbol==inv_gene & wynik$enhancer_pvalue<0.01 & !is.na(wynik$enhancer_pvalue)]
  enh_stop <- enh_stop[1]
  enh_data_to_calculate <- enhancery[enhancery$start==enh_start & enhancery$end==enh_stop,]
  x <- c(x,as.numeric(data[data$symbol==inv_gene,patients]))
  y <- c(y, as.numeric(enh_data_to_calculate[,patients]))
  kor <- cor.test(as.numeric(data[data$symbol==inv_gene,patients]), as.numeric(enh_data_to_calculate[,patients]), method='spearman')
  kor_R <- c(kor_R, kor$estimate)
  kor_p <- c(kor_p, kor$p.value)
}
median(kor_R) # 0.45
median(kor_p) # 0.04
sd(kor_R)

# 
# # Only PA:
# patients_sub <- patients[grepl('PA', patients)]
# kor_R_1 <- kor_p_1 <-c()
# x <- y <- c()
# for (gene in unique(geny_kontakty_rozn$Var1)) {
#   inv_gene <- gene
#   enh_start <- wynik$enhancer_start[wynik$symbol==inv_gene & wynik$enhancer_pvalue<0.01 & !is.na(wynik$enhancer_pvalue)]
#   enh_start <- enh_start[1]
#   enh_stop <- wynik$enhancer_stop[wynik$symbol==inv_gene & wynik$enhancer_pvalue<0.01 & !is.na(wynik$enhancer_pvalue)]
#   enh_stop <- enh_stop[1]
#   #enhancery[enhancery$start==enh_start & enhancery$end==enh_stop+1,]
#   enh_data_to_calculate <- enhancery[enhancery$start==enh_start & enhancery$end==enh_stop,]
#   x <- c(x,as.numeric(data[data$symbol==inv_gene,patients_sub]))
#   y <- c(y, as.numeric(enh_data_to_calculate[,patients_sub]))
#   kor <- cor.test(as.numeric(data[data$symbol==inv_gene,patients_sub]), as.numeric(enh_data_to_calculate[,patients_sub]), method='spearman')
#   kor_R_1 <- c(kor_R_1, kor$estimate)
#   kor_p_1 <- c(kor_p_1, kor$p.value)
# }
# 
# median(na.omit(kor_R_1)) # 0.4
# median(na.omit(kor_p_1)) # 0.42
# sd(kor_R_1, na.rm = T)
# 
# # Only DA:
# patients_sub <- patients[grepl('DA', patients)]
# kor_R_23 <- kor_p_23 <-c()
# x <- y <- c()
# for (gene in unique(geny_kontakty_rozn$Var1)) {
#   inv_gene <- gene
#   enh_start <- wynik$enhancer_start[wynik$symbol==inv_gene & wynik$enhancer_pvalue<0.01 & !is.na(wynik$enhancer_pvalue)]
#   enh_start <- enh_start[1]
#   enh_stop <- wynik$enhancer_stop[wynik$symbol==inv_gene & wynik$enhancer_pvalue<0.01 & !is.na(wynik$enhancer_pvalue)]
#   enh_stop <- enh_stop[1]
#   #enhancery[enhancery$start==enh_start & enhancery$end==enh_stop+1,]
#   enh_data_to_calculate <- enhancery[enhancery$start==enh_start & enhancery$end==enh_stop,]
#   x <- c(x,as.numeric(data[data$symbol==inv_gene,patients_sub]))
#   y <- c(y, as.numeric(enh_data_to_calculate[,patients_sub]))
#   kor <- cor.test(as.numeric(data[data$symbol==inv_gene,patients_sub]), as.numeric(enh_data_to_calculate[,patients_sub]), method='spearman')
#   kor_R_23 <- c(kor_R_23, kor$estimate)
#   kor_p_23 <- c(kor_p_23, kor$p.value)
# }
# 
# median(na.omit(kor_R_23)) # 0.43
# median(na.omit(kor_p_23)) # 0.27
# sd(kor_R_23, na.rm = T) #0.45
# 
# # Only GB|PG:
# patients_sub <- patients[grepl('GB|PG', patients)]
# kor_R_4 <- kor_p_4 <-c()
# x <- y <- c()
# for (gene in unique(geny_kontakty_rozn$Var1)) {
#   inv_gene <- gene
#   enh_start <- wynik$enhancer_start[wynik$symbol==inv_gene & wynik$enhancer_pvalue<0.01 & !is.na(wynik$enhancer_pvalue)]
#   enh_start <- enh_start[1]
#   enh_stop <- wynik$enhancer_stop[wynik$symbol==inv_gene & wynik$enhancer_pvalue<0.01 & !is.na(wynik$enhancer_pvalue)]
#   enh_stop <- enh_stop[1]
#   #enhancery[enhancery$start==enh_start & enhancery$end==enh_stop+1,]
#   enh_data_to_calculate <- enhancery[enhancery$start==enh_start & enhancery$end==enh_stop,]
#   x <- c(x,as.numeric(data[data$symbol==inv_gene,patients_sub]))
#   y <- c(y, as.numeric(enh_data_to_calculate[,patients_sub]))
#   kor <- cor.test(as.numeric(data[data$symbol==inv_gene,patients_sub]), as.numeric(enh_data_to_calculate[,patients_sub]), method='spearman')
#   kor_R_4 <- c(kor_R_4, kor$estimate)
#   kor_p_4 <- c(kor_p_4, kor$p.value)
# }
# median(na.omit(kor_R_4)) # 0.23
# median(na.omit(kor_p_4)) # 0.38
# sd(kor_R_4, na.rm = T)

# Plotting acetylations vs expression quartiles:

data$grade1_mean <- apply(data[,grepl('PA', colnames(data))], 1, mean)
data$grade23_mean <- apply(data[,grepl('DA', colnames(data))], 1, mean)
data$grade4_mean <- apply(data[,grepl('GB|PG', colnames(data))], 1, mean)
data$mean_all <- apply(data[,grepl('PA|DA|GB|PG', colnames(data))], 1, mean)
data$grade1_expr_quartile <- as.factor(ntile(data$grade1_mean, 4)) #Dodaje info o kwartylu ekspresji
data$grade23_expr_quartile <- as.factor(ntile(data$grade23_mean, 4))
data$grade4_expr_quartile <- as.factor(ntile(data$grade4_mean, 4))

wynik_istotne <- wynik[wynik$enhancer_pvalue<0.01 & !is.na(wynik$enhancer_pvalue),]

for (gene in data$symbol) {
  chr <- unique(wynik_istotne$chr[wynik_istotne$symbol==gene])
  enh_start <- unique(wynik_istotne$enhancer_start[wynik_istotne$symbol==gene])[1] #If there are more enhancers I took the first one (it should be most significant)
  data$grade1_acetylation[data$symbol==gene] <-
    apply(enhancery[enhancery$chr==chr & enhancery$start==enh_start,grepl('PA', colnames(enhancery))], 1, mean)
  data$grade23_acetylation[data$symbol==gene] <-
    apply(enhancery[enhancery$chr==chr & enhancery$start==enh_start,grepl('DA', colnames(enhancery))], 1, mean)
  data$grade4_acetylation[data$symbol==gene] <-
    apply(enhancery[enhancery$chr==chr & enhancery$start==enh_start,grepl('GB|PGA', colnames(enhancery))], 1, mean)
}

cor.test(c(data$grade1_mean, data$grade23_mean, data$grade4_mean), c(data$grade1_acetylation, data$grade23_acetylation, data$grade4_acetylation))
gr1_cor <- cor.test(data$grade1_mean, data$grade1_acetylation, method='spearman') # grade I
gr23_cor <- cor.test(data$grade23_mean, data$grade23_acetylation, method='spearman') # grade23
gr4_cor <- cor.test(data$grade4_mean, data$grade4_acetylation, method='spearman') #grade4
gr1_cor$estimate
gr1_cor$p.value
gr23_cor$estimate
gr23_cor$p.value
gr4_cor$estimate
gr4_cor$p.value

ggplot(data, aes(x=grade1_expr_quartile, y=grade1_acetylation)) +
  geom_boxplot() +
  xlab('Grade I \nexpression quartiles') + ylab ('Grade I enhancer H3K27ac levels' ) + ylim(c(0,15)) +
  annotate("text", x = 3, y = 14, label = paste0('rho = ', round(gr1_cor$estimate,2))) +
  geom_point(color='#56B4E9', size=4, alpha=0.6) + #size=2, color = 'dimgrey'
  theme(text = element_text(size=15))
ggsave('Grade1_expression_quartiles_vs_H3K27ac.png', width=2.5, height = 5)

ggplot(data, aes(x=grade23_expr_quartile, y=grade23_acetylation)) +
  geom_boxplot() +
  xlab('Grade II/III \nexpression quartiles') + ylab ('Grade II/III enhancer H3K27ac levels' ) + ylim(c(0,15)) +
  annotate("text", x = 3, y = 14, label = paste0('rho = ', round(gr23_cor$estimate, 2))) +
  geom_point(color = '#E69F00', size=4, alpha=0.6) +
  theme(text = element_text(size=15))
ggsave('Grade23_expression_quartiles_vs_H3K27ac.png', width=2.5, height = 5)

ggplot(data, aes(x=grade4_expr_quartile, y=grade4_acetylation)) +
  geom_boxplot() +
  xlab('Grade IV \nexpression quartiles') + ylab ('Grade IV enhancer H3K27ac levels' ) + ylim(c(0,15)) +
  annotate("text", x = 3, y = 14, label = paste0('rho = ', round(gr4_cor$estimate, 2))) +
  geom_point(color = '#999999', size=4, alpha=0.6) +
  theme(text = element_text(size=15))
ggsave('Grade4_expression_quartiles_vs_H3K27ac.png', width=2.5, height = 5)

### PERMUTATION

gr1_cor_median <- gr23_cor_median <- gr4_cor_median <- c()
for (rep in 1:100) {
  random_rows<- sample(1:nrow(data), nrow(data))
gr1_cor <- cor.test(data$grade1_mean, data$grade1_acetylation[random_rows], method='spearman') # grade I
gr23_cor <- cor.test(data$grade23_mean, data$grade23_acetylation[random_rows], method='spearman') # grade23
gr4_cor <- cor.test(data$grade4_mean, data$grade4_acetylation[random_rows], method='spearman') #grade4
gr1_cor_median <- c(gr1_cor_median, gr1_cor$estimate)
gr23_cor_median <- c(gr23_cor_median, gr23_cor$estimate)
gr4_cor_median <- c(gr4_cor_median, gr4_cor$estimate)
}
max(gr1_cor_median)
min(gr1_cor_median)
median(gr1_cor_median)
sd(gr1_cor_median)

max(gr23_cor_median)
min(gr23_cor_median)
median(gr23_cor_median)
sd(gr23_cor_median)

max(gr4_cor_median)
min(gr4_cor_median)
median(gr4_cor_median)
sd(gr4_cor_median)

##### Enhancer methylations (DNA methylation):
system("bash ./long_range_interactions/enhancer_methylations/1_run-python_enh_methylation_means.sh")
setwd('./long_range_interactions/enhancer_methylations/')
system("python ./2_skladanie_pacjentow_w_calosc.py")
setwd('../..')
getwd()
# Read-in the enhancer methylation data and find differentially methylated ones:
enhancer_methylations <- read.csv('./long_range_interactions/enhancer_methylations/enhancers_methylations.txt', header=T, sep='\t', stringsAsFactors = F)
for (ln in 1:nrow(enhancer_methylations)) {
  if (sum(!is.na(enhancer_methylations[ln, grepl('PA', colnames(enhancer_methylations))]))>0 & sum(!is.na(enhancer_methylations[ln, grepl('GB|PG', colnames(enhancer_methylations))]))   ) {
    test <- (wilcox.test(as.numeric(enhancer_methylations[ln, grepl('PA', colnames(enhancer_methylations))]), as.numeric(enhancer_methylations[ln, grepl('GB|PG', colnames(enhancer_methylations))])))
    enhancer_methylations$wilcox_p[ln] <- test$p.value
  }
}

# Finding genes contacting enhancers being sign. methylated (DNA):
enh_meth_sign <- enhancer_methylations[enhancer_methylations$wilcox_p<0.01 & !is.na(enhancer_methylations$wilcox_p),]
wynik <- wynik[!is.na(wynik$enhancer_start),]
wynik$wilcox_p_methyl <- NA
for (i in 1:nrow(wynik)) {
  for (j in 1:nrow(enh_meth_sign)){
    if (wynik$chr[i] == enh_meth_sign$chr[j] & wynik$enhancer_start[i] == enh_meth_sign$enh_start[j] ){ 
      wynik$wilcox_p_methyl[i] <-enh_meth_sign$wilcox_p[j]
    }
  }
}
# How many genes contact diff. meth enhancers?
length(unique(wynik$symbol[!is.na(wynik$wilcox_p_methyl) & wynik$wilcox_p_methyl<0.01]))

# Adding to the object containing genes with significant enhancer acteylations significant methylations:
enh_meth_sign <- enhancer_methylations[enhancer_methylations$wilcox_p<0.01 & !is.na(enhancer_methylations$wilcox_p),]
wynik_istotne<- wynik[wynik$enhancer_pvalue<0.01&wynik$wilcox_p_methyl<0.01,]
wynik_istotne<- na.omit(wynik_istotne)

# Finding genes contacting enhancers being both significantly acetylated (histone) and sign. methylated (DNA):
genes_sign_meth_acet <- wynik_istotne[wynik_istotne$wilcox_p_methyl<0.01 & !is.na(wynik_istotne$wilcox_p_methyl),]
length(unique(genes_sign_meth_acet$DE_ENSEMBL))
write.csv(genes_sign_meth_acet, 'genes_sign_meth_acet.csv')
#write.csv(genes_sign_meth_acet, 'genes_sign_meth_acet.rnk') ########## robic to?

# Calculate correlation of methylation and acetylation
# Remove duplicate rows of the dataframe using symbol variable:
patients <- colnames(enhancery)[grepl('PA|GB|PG',colnames(enhancery))] #enh_meth_sign]
patients <- patients [patients %in% colnames(enh_meth_sign)]

genes_sign_meth_acet <- distinct(genes_sign_meth_acet,symbol, .keep_all= TRUE)
plot_acet_meth <- genes_sign_meth_acet[,c('DE_ENSEMBL', 'symbol', 'chr', 'enhancer_start', 'enhancer_stop')]
plot_acet_meth <- cbind(plot_acet_meth, data.frame(matrix(ncol = length(patients), nrow = nrow(plot_acet_meth))))
colnames(plot_acet_meth) <- c('DE_ENSEMBL', 'symbol', 'chr', 'enhancer_start', 'enhancer_stop', patients)
tmp_acet <- plot_acet_meth

korelacje_R <- korelacje_p <- c()
for (line in 1:nrow(genes_sign_meth_acet)) {
  chr <- genes_sign_meth_acet$chr[line]
  enh_start <- genes_sign_meth_acet$enhancer_start[line]
  korelacja <- cor.test(
    as.numeric(enhancery[enhancery$chr==chr & enhancery$start==enh_start, patients]) , 
    as.numeric(enh_meth_sign[enh_meth_sign$chr==chr & enh_meth_sign$enh_start==enh_start, patients]) ,
    method='spearman'
  )
    korelacje_R <- c(korelacje_R, korelacja$estimate)
    korelacje_p <- c(korelacje_p, korelacja$p.value)
    genes_sign_meth_acet$kor_R[line] <- korelacja$estimate
    genes_sign_meth_acet$kor_p[line] <- korelacja$p.value
  plot_acet_meth [plot_acet_meth$chr == chr & plot_acet_meth$enhancer_start == enh_start , patients] <- enh_meth_sign[enh_meth_sign$chr==chr & enh_meth_sign$enh_start==enh_start, patients]
  tmp_acet [plot_acet_meth$chr == chr & plot_acet_meth$enhancer_start == enh_start , patients] <- enhancery[enhancery$chr==chr & enhancery$start==enh_start, patients]
}
write.csv(genes_sign_meth_acet, 'genes_sign_meth_acet.csv') # suppl Table
median(korelacje_R)
sd(korelacje_R)
median(korelacje_p)
plot_acet_meth <- reshape2::melt(plot_acet_meth, id.vars=c("DE_ENSEMBL", "symbol", 'chr', 'enhancer_start', 'enhancer_stop'),
                        measure.vars=c(patients),
                        variable.name = "samples",
                        value.name="enhancer_DNA_methylation")
tmp_acet <- reshape2::melt(tmp_acet, id.vars=c("DE_ENSEMBL", "symbol", 'chr', 'enhancer_start', 'enhancer_stop'),
                       measure.vars=c(patients),
                       variable.name="samples",
                       value.name="enhancer_H3K27acetylation")
plot_acet_meth$enhancer_H3K27acetylation <- tmp_acet$enhancer_H3K27acetylation
plot_acet_meth$Grade[grepl('PA',plot_acet_meth$samples)]<-'grade_I'
plot_acet_meth$Grade[grepl('DA',plot_acet_meth$samples)]<-'grade_II/III'
plot_acet_meth$Grade[grepl('GB|PG',plot_acet_meth$samples)]<-'grade_IV'
cor.test(plot_acet_meth$enhancer_DNA_methylation , plot_acet_meth$enhancer_H3K27acetylation, method='spearman') # rho = -0.26 p= 0.008

ggplot(plot_acet_meth, aes(x=enhancer_DNA_methylation, y=enhancer_H3K27acetylation, color=Grade, shape=symbol)) +
geom_point( size=4, alpha=0.8) +
  xlab('Enhancer DNA CpG average methylation [%]') +
  ylab('Enhancer H3K27ac coverage') +
  guides(color=guide_legend(title="WHO glioma grade")) +
  guides(shape=guide_legend(title="Gene symbol")) +
  theme(text = element_text(size=15)) +
  scale_color_manual(values=c('#56B4E9', '#999999')) +
  scale_shape_manual(values=c(15, 2, 17, 18, 19, 1, 7, 8, 6)) +
  scale_y_sqrt() +
  scale_x_sqrt()
ggsave('Enhancer_H3K27ac_DNAmethylation.png', width = 7, height = 5)

# Corrplots of correlations between these 9 genes' expressions and H3K27acetylation or DNA methylation:
# Corrplot for all PCDHG genes expression and contacting enhancers acetylation:
library('corrplot')
library(magrittr)
# 1) Expression and contacting enhancers acetylation
genes_9 <- genes_sign_meth_acet$symbol
genes_9_expressions <- data.frame(t(data[data$symbol %in% genes_9, c('symbol', pacjenci)]))
colnames(genes_9_expressions) <- genes_9_expressions[rownames(genes_9_expressions)=='symbol',]
genes_9_expressions <- genes_9_expressions[2:nrow(genes_9_expressions),] 
genes_9_expressions[,] %<>% lapply(function(x) as.numeric(as.character(x)))

genes_9_string <- 'EGFR|GLI3|ITGB3BP|FARSB|DNAH11|ADGRA1|PHACTR1|EDN1|ADGRF5P1'
genes_9_enhancers_region_start <- unique(wynik_istotne$region_start[!is.na(wynik_istotne$enhancer_start) & grepl(genes_9_string, wynik_istotne$symbol)])
genes_9_enhancers_region_stop <- unique(wynik_istotne$region_stop[!is.na(wynik_istotne$enhancer_stop) & grepl(genes_9_string, wynik_istotne$symbol)])
genes_9_enhancers_region_chr <- unique(wynik_istotne$chr[!is.na(wynik_istotne$enhancer_start) & grepl(genes_9_string, wynik_istotne$symbol)])
genes_9_enhancers <- data.frame(t(enhancery[enhancery$chr %in% genes_9_enhancers_region_chr & 
                                                  enhancery$start %in% genes_9_enhancers_region_start, c( 'start', pacjenci)]))
colnames(genes_9_enhancers) <- paste0(genes_9_enhancers_region_chr, '_', 
                                      genes_9_enhancers_region_start, '_', 
                                      genes_9_enhancers_region_stop)
genes_9_enhancers <- genes_9_enhancers[2:nrow(genes_9_enhancers),] 
genes_9_enhancers[,] %<>% lapply(function(x) as.numeric(as.character(x)))

corr_mat=cor(genes_9_expressions, genes_9_enhancers, method="s")
corrplot(corr_mat,tl.col = 'black', col = cor_col(20))

png('corplot_for_9_genes_and_9_enhancers_acetylation_expression.png')
corrplot(cor_for_plot,tl.col = 'black', col = cor_col(20))
dev.off()

# 2) Expression and contacting enhancers DNA methylation
pacjenci <- colnames(enhancer_methylations)[grepl('PA|GB|PG',colnames(enhancer_methylations))] #enh_meth_sign]
genes_9 <- genes_sign_meth_acet$symbol
genes_9_expressions <- data.frame(t(data[data$symbol %in% genes_9, c('symbol', pacjenci)]))
colnames(genes_9_expressions) <- genes_9_expressions[rownames(genes_9_expressions)=='symbol',]
genes_9_expressions <- genes_9_expressions[2:nrow(genes_9_expressions),] 
genes_9_expressions[,] %<>% lapply(function(x) as.numeric(as.character(x)))

genes_9_string <- 'EGFR|GLI3|ITGB3BP|FARSB|DNAH11|ADGRA1|PHACTR1|EDN1|ADGRF5P1'
genes_9_enhancers_region_start <- unique(wynik_istotne$region_start[!is.na(wynik_istotne$enhancer_start) & grepl(genes_9_string, wynik_istotne$symbol)])
genes_9_enhancers_region_stop <- unique(wynik_istotne$region_stop[!is.na(wynik_istotne$enhancer_stop) & grepl(genes_9_string, wynik_istotne$symbol)])
genes_9_enhancers_region_chr <- unique(wynik_istotne$chr[!is.na(wynik_istotne$enhancer_start) & grepl(genes_9_string, wynik_istotne$symbol)])
genes_9_enhancers <- data.frame(t(enhancer_methylations[enhancer_methylations$chr %in% genes_9_enhancers_region_chr & 
                                      enhancer_methylations$enh_start %in% genes_9_enhancers_region_start, c('enh_start', pacjenci)]))
colnames(genes_9_enhancers) <- paste0(genes_9_enhancers_region_chr, '_', 
                                      genes_9_enhancers_region_start, '_', 
                                      genes_9_enhancers_region_stop)
genes_9_enhancers <- genes_9_enhancers[2:nrow(genes_9_enhancers),] 
genes_9_enhancers[,] %<>% lapply(function(x) as.numeric(as.character(x)))
corr_mat=cor(genes_9_expressions, genes_9_enhancers, method="spearman")
corrplot(corr_mat)

corr_mat2=cor(genes_9_expressions, genes_9_enhancers, method="pearson")
corrplot(corr_mat2)
png('corplot_for_9_genes_and_9_enhancers_DNAmethylation_expression.png')
corrplot(cor_for_plot)
dev.off()