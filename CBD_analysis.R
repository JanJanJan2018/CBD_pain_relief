# CBD GEO data sets inquiry

GPL6244 <- read.delim('GPL6244-17930.txt', sep='\t', quote="", header=TRUE,
                  comment.char='#', na.strings=c('','NA','---'))

GPL4133 <- read.delim('GPL4133-12599.txt', sep='\t', quote="", header=TRUE,
                       comment.char='#', na.strings=c('','NA','---'))

GSE57571 <- read.delim('GSE57571_series_matrix.txt',sep='\t', quote="", header=TRUE,
                        comment.char='!', na.strings=c('','NA'))
colnames(GSE57571) <- gsub('X.', '',colnames(GSE57571))

GSE57978 <- read.delim('GSE57978_series_matrix.txt', sep='\t', quote="", header=TRUE,
                        comment.char='!', na.strings=c('','NA'))
colnames(GSE57978) <- gsub('X.', '', as.character(colnames(GSE57978)))

GSE57978merged <- merge(GSE57978,GPL6244, by.x='ID_REF.', by.y='ID')
colnames(GSE57978merged)[c(2,4,6)] <- paste(as.character(colnames(GSE57978merged)[c(2,4,6)]), 
                                                         'CBD', sep='')
colnames(GSE57978merged)[c(3,5,7)] <- paste(as.character(colnames(GSE57978merged)[c(3,5,7)]), 
                                            'CTRL', sep='')

Sym <- strsplit(as.character(GSE57978merged$gene_assignment), '//')
Ref_Seq <- lapply(Sym,'[',1)
Symbol <- lapply(Sym,'[',2)
HGNC <- lapply(Sym,'[',3)

Ref_Seq <- t(data.frame(Ref_Seq))
row.names(Ref_Seq) <- NULL
colnames(Ref_Seq) <- 'Ref_Seq'

Symbol <- t(data.frame(Symbol))
row.names(Symbol) <- NULL
colnames(Symbol) <- 'Symbol'

HGNC <- t(data.frame(HGNC))
row.names(HGNC) <- NULL
colnames(HGNC) <- 'HGNC'

names <- cbind(Ref_Seq,Symbol,HGNC)

GSE56978platform <- cbind(GSE57978merged, names)
gse56978 <- GSE56978platform[complete.cases(GSE56978platform),]




GSE57571merged <- merge(GSE57571, GPL4133, by.x='ID_REF.', by.y='ID')
GSE57571merged <- GSE57571merged[,-c(26)]
gse57571 <- GSE57571merged[complete.cases(GSE57571merged),]
colnames(gse57571)[5] <- paste(colnames(gse57571)[5], 'CTRL', sep='')
colnames(gse57571)[3:4] <- paste(colnames(gse57571)[3:4],'REP2', sep='')
colnames(gse57571)[6] <- paste(colnames(gse57571)[6], 'CBD', sep='')
colnames(gse57571)[2] <- paste(colnames(gse57571)[2], 'REP1', sep='')
colnames(gse57571)[1] <- 'ID_REF'

write.csv(gse56978, 'gse56978_full_set.csv', row.names=FALSE)
write.csv(gse57571, 'gse57571_full_set.csv', row.names=FALSE)

###############################################################################

gse56978 <- read.csv('gse56978_full_set.csv',  header=TRUE, sep=',')

#divide this large data set into four parts 22152/4=5538 rows long each
# 1:5538, 5539:11076, 11077:16614, 16615:22152
gse56978_1 <- gse56978[1:5538,]
gse56978_2 <- gse56978[5539:11076,]
gse56978_3 <- gse56978[11077:16614,]
gse56978_4 <- gse56978[16615:21152,]

write.csv(gse56978_1, 'gse56978_1.csv', row.names=FALSE)
write.csv(gse56978_2, 'gse56978_2.csv', row.names=FALSE)
write.csv(gse56978_3, 'gse56978_3.csv', row.names=FALSE)
write.csv(gse56978_4, 'gse56978_4.csv', row.names=FALSE)

gse57571 <- read.csv('gse57571_full_set.csv',  header=TRUE, sep=',')

gliomaCBD <- gse56978[,c(2:7,20)]
sebumCBD <- gse57571[c(2:6,15)]

gliomaCBD <- gliomaCBD[!duplicated(gliomaCBD$Symbol),]
sebumCBD <- sebumCBD[!duplicated(sebumCBD$GENE_SYMBOL),]

namesS <- sebumCBD$GENE_SYMBOL
namesG <- gliomaCBD$Symbol

row.names(sebumCBD) <- namesS
row.names(gliomaCBD) <- namesG

gliomaCBD <- gliomaCBD[,-7]
sebumCBD <- sebumCBD[,-6]

write.csv(gliomaCBD, 'gliomaCBD.csv', row.names=TRUE)
write.csv(sebumCBD, 'sebumCBD.csv', row.names=TRUE)

# note that the sebum samples show the difference between the control vehicle and 
# CBD applied 24 hours in rep2 and rep3, thus the negative values in samples GSM1385014-GSM1385017
# The two data sets cannot be combined due to this fact, even though both are array data
# all of these samples are CBD treated as the control variable is FALSE for all observations

# Gene Symbols referenced and tested in article on GSE57571 sebacious:
# TNFA,  TRPV1, TRPV2, TRPV4, IL1B, IL6, NRIP1, TRIB3 ( SINK), ARHGAP9, ZM,  MKI67 

#none of these commented are in this data
#TNFA <- grep('TNFA', namesS)
#TRPV1 <- grep('TRPV1', namesS)
#TRPV4 <- grep('TRPV4', namesS)
#ARHGAP9 <- grep('ARHGAP9', namesS)
#ZM <- grep('ZM', namesS)

MKI67 <- grep('MKI67', namesS)
MKI67 <- MKI67[2]
TRPV2 <- grep('TRPV2', namesS)
IL1B <- grep('IL1B', namesS)
IL6 <- grep('IL6', namesS)
IL6 <- IL6[2]
NRIP1 <- grep('NRIP1', namesS)
NRIP1 <- NRIP1[1]
TRIB3 <- grep('TRIB3', namesS)

mki67 <- sebumCBD[MKI67,]
trpv2 <- sebumCBD[TRPV2,]
il1b <- sebumCBD[IL1B,]
il6 <- sebumCBD[IL6,]
nrip1 <- sebumCBD[NRIP1,]
trib3 <- sebumCBD[TRIB3,]
CNR1s <- sebumCBD[grep('CNR1', namesS),]
CNR1s <- CNR1s[1,]
CNRIP1s <- sebumCBD[grep('CNRIP1', namesS),]
CNR2s <- sebumCBD[grep('CNR2', namesS),]

genesSebumCBD <- rbind(mki67, trpv2, il1b, il6, nrip1, trib3,CNR1s,CNRIP1s,CNR2s)
write.csv(genesSebumCBD,'genesSebumCBD.csv', row.names=TRUE)

# all paraphrased summaries from gene summaris on ncbi.nlm.nih.gov
TNFA <- gliomaCBD[grep('TNFA', namesG),] #a cytokine that binds to specific TNF receptors,
#involved in cell proliferation, cell apoptosis, cell differentiation, lipid metabolism, and 
#coagulation. It is associated with cancer, autoimmune diseases, and diabetes, in mice it has
#a neuroprotective property
S100B <- gliomaCBD[grep('S100B', namesG),]# a protein involved with 13 members that are shown
# to be involved in cell cycle progression and differentiation, and mutations can be involved in 
# neurological and neoplastic diseases like AD, epilepsy, Type I diabetes, melanoma, and ALS
IL1B <- gliomaCBD[grep('IL1B', namesG),]#associates with pain and inflammation mediator and
#processed by CASP1 caspase 1 to full active form, IL1B is a cytokine, responds to macrophages,
#and is found in CNS for pain hypersensitivity when cyclooxygenase-2 (PTGS2/COX2) is present, IL1B
#also is involved in cell proliferation, cell apoptosis, and differentiation
IL6 <- gliomaCBD[grep('IL6', namesG),] #associated with local inflammation and cytokine production

CNRIP1 <- gliomaCBD[grep('CNRIP1', namesG),]
CNR1 <- gliomaCBD[grep('CNR1', namesG),]
CNR1 <- CNR1[2,]#protein coding, cannabinoid receptor 1 gene
CNR2 <- gliomaCBD[grep('CNR2', namesG),]#cannabinoid receptor 2 gene, 
#involved in CNS effects by psychoactive properties of cannabinoids for mood and cognition


#combine these genes to a table of inflammation and disease genes for gliomaCBD

diseaseGenesGlioma <- rbind(TNFA,S100B,IL1B,IL6,CNRIP1,CNR1,CNR2)
write.csv(diseaseGenesGlioma, 'GliomaDiseaseGenes.csv', row.names=TRUE)

##############################################################################################

gliomaCBD <- read.csv('gliomaCBD.csv', row.names=1, header=TRUE, sep=',')
sebumCBD <- read.csv('sebumCBD.csv', row.names=1, header=TRUE, sep=',')

namesG <- row.names(gliomaCBD)
namesG <- gsub(' ','', namesG)
gliomaCBD$Symbol <- namesG

namesS <- row.names(sebumCBD)
sebumCBD$Symbol <- namesS

# I made a summary of the genes interested in with inflammation, cell prolif/death, 
# tumor suppr/inducers saved as 'NCBI_Ref_Seq_Gene_SUMM.csv' to add to CBD data set

summ <- read.delim('NCBI_Ref_Seq_Gene_SUMM.csv', sep=',', comment.char='#', skip =1)

# merge to gene summaries
gliomaCBD_summ <- merge(summ, gliomaCBD, by.x='Gene', by.y='Symbol')
sebumCBD_summ <- merge(summ, sebumCBD, by.x='Gene', by.y='Symbol')

write.csv(gliomaCBD_summ, 'gliomaCBD_summ.csv', row.names=FALSE)
write.csv(sebumCBD_summ, 'sebumCBD_summ.csv', row.names=FALSE)

library(dplyr)

# add a DE column and fold change column to the sebumCBD_summ data
DE_sebum <- mutate(sebumCBD_summ, DE_CBD_mns_CTRL = GSM1385017.CBD-GSM1385016.CTRL)
mag_sebum <- mutate(DE_sebum, magnitude_CBD = abs(GSM1385017.CBD-GSM1385016.CTRL))
FC_sebum <- mutate(mag_sebum, Fold_Change_CBD_2_CTRL=GSM1385017.CBD/GSM1385016.CTRL)

FC_sebum_DE_order <- FC_sebum[order(FC_sebum$DE_CBD_mns_CTRL,decreasing=TRUE),]
FC_sebum_FC_order <- FC_sebum[order(FC_sebum$Fold_Change_CBD_2_CTRL, decreasing=TRUE),]

write.csv(FC_sebum_FC_order, 'FC_sebum_order.csv', row.names=FALSE)
write.csv(FC_sebum_DE_order, 'DE_sebum_order.csv', row.names=FALSE)

FC_sebum_order <- read.csv('FC_sebum_order.csv', sep=',', header=TRUE)

library(dplyr)

DE <- FC_sebum_order[c(1:4,20:22),]
DE2 <- mutate(DE, expression=DE$DE_CBD_mns_CTRL>0)
DE2$expression <- gsub('TRUE','up-regulated',DE2$expression)
DE2$expression <- gsub('FALSE','down-regulated', DE2$expression)

library(ggplot2)

png('GenesRegulatedMostSebumCBD.png', width=768, height=576)
g <- ggplot(DE2, aes(x=as.factor(Gene), y=Fold_Change_CBD_2_CTRL))
g= g+xlab('Genes in CBD Treated Sebum to Non-Treated')
g= g+ ylab('Fold Change CBD Treated Sebum to Non-Treated')
g= g+ geom_point(aes(colour=expression),size=6, alpha=0.9)
g
dev.off()

###################################################################################################3

FC_sebum_FC_order <- read.csv('FC_sebum_order.csv', sep=',', header=TRUE)

# pull in the data set on all genes from the UL target risk gene study of 121 samples
# 51 non and 70 UL, 'all_common_12173_130_fold_magnitude.csv' and get the genes used in
# the CBD data on inflammation, cytokines, pain, hormones, etc.

UL_FC <- read.csv('all_common_12173_130_fold_magnitude.csv', sep=',', header=TRUE, row.names=1)

FC_cbd <- data.frame(FC_sebum_FC_order[,1:3])
write.csv(FC_cbd,'hormon-pain-genes.csv', row.names=FALSE)

UL_inf_pain_hormone <- merge(FC_cbd,UL_FC, by.x='Gene', by.y='GENE_SYMBOL')

write.csv(UL_inf_pain_hormone, 'UL_pain_hormone_genes.csv', row.names=FALSE)

# compare those genes expressed in UL compared to CBD treated sebum array samples for fold change
# and differential expression

###################################################################################################3

# now compare the CBD and CTRL groups of the glioma CBD treated samples, 3 each

CBD_glioma <- read.csv('gliomaCBD_summ.csv', sep=',', header=TRUE)

CBD_trtd <- CBD_glioma[,c(4,6,8)]
names <- as.character(CBD_glioma$Gene)
row.names(CBD_trtd) <- names

CBD_cntrl <- CBD_glioma[,c(5,7,9)]
row.names(CBD_cntrl) <- names

gliomaMeans <- data.frame(rowMeans(CBD_cntrl))
colnames(gliomaMeans) <- 'Control_Mean'

gliomaCBDMeans <- data.frame(rowMeans(CBD_trtd))
colnames(gliomaCBDMeans) <- 'CBD_Treated_Mean'

gliomaCBD_stats <- cbind(gliomaMeans, gliomaCBDMeans)

library(dplyr)

gliomaStats <- mutate(gliomaCBD_stats, DifferenceInMeansFromCBD=CBD_Treated_Mean-Control_Mean)
gliomaStats2 <- mutate(gliomaStats, FoldChange=CBD_Treated_Mean/Control_Mean)

CBD_glioma_stats <- cbind(CBD_glioma, gliomaStats2)

write.csv(CBD_glioma_stats, 'CBD_glioma_stats.csv', row.names=FALSE)
CBD_glioma_stats <- read.csv('CBD_glioma_stats.csv', sep=',', header=TRUE)

CBD_glioma_FC_ordered <- CBD_glioma_stats[order(CBD_glioma_stats$FoldChange,
                                                decreasing=TRUE),]
write.csv(CBD_glioma_FC_ordered,'CBD_glioma_FC_ordered.csv', row.names=FALSE)


CBD_glioma_DE_ordered <- CBD_glioma_stats[order(CBD_glioma_stats$DifferenceInMeansFromCBD,
                                                decreasing=TRUE),]

write.csv(CBD_glioma_DE_ordered,'CBD_glioma_DE_ordered.csv', row.names=FALSE)

library(dplyr)

DE <- CBD_glioma_FC_ordered[c(1:2,27:30),]
DE2 <- mutate(DE, expression=DE$DifferenceInMeansFromCBD>0)
DE2$expression <- gsub('TRUE','up-regulated',DE2$expression)
DE2$expression <- gsub('FALSE','down-regulated', DE2$expression)

library(ggplot2)

png('GenesRegulatedMostGliomaCBD.png', width=768, height=576)
g <- ggplot(DE2, aes(x=as.factor(Gene), y=FoldChange))
g= g+xlab('Genes in CBD Treated Glioma to Non-Treated')
g= g+ ylab('Fold Change CBD Treated Glioma to Non-Treated')
g= g+ geom_point(aes(colour=expression),size=6, alpha=0.9)
g
dev.off()


###################################################################################################3

# set the directory for each of these files then go to the main directory

painGenes <- read.csv('GenesSumm1.csv', header=TRUE, sep=',')
colnames(painGenes)[1] <- 'Gene'

EndocrineGenes <- read.csv('GenesSumm2.csv', header=TRUE, sep=',')
colnames(EndocrineGenes)[1] <- 'Gene'

# go to main directory of this script and files

# retrieve the files from the folder of the pancreatic cancer to compare genes
GSE131859 <- read.delim('GSE131859_non-normalized.txt',sep='\t', quote="", header=TRUE,
                        comment.char='!', na.strings=c('','NA'))

GPL14951 <- read.delim('GPL14951-11332.txt', sep='\t', quote="", header=TRUE,
                      comment.char='#', na.strings=c('','NA','---'))
#this file above is 60 mb and broken into smaller sections to combine from github in 3 parts a:c

GPL14951_a <- GPL14951[1:9792,]
write.csv(GPL14951_a, 'GPL14951_a.csv', row.names=FALSE)
GPL14951_b <- GPL14951[9793:19584,]
write.csv(GPL14951_b, 'GPL14951_b.csv', row.names=FALSE)
GPL14951_c <- GPL14951[19585:29377,]
write.csv(GPL14951_c, 'GPL14951_c.csv', row.names=FALSE)

GPL14951_a <- read.csv('GPL14951_a.csv', sep=',',  header=TRUE,
                        na.strings=c('','NA','---'))
GPL14951_b <- read.csv('GPL14951_b.csv', sep=',', header=TRUE,
                          na.strings=c('','NA','---'))
GPL14951_c <- read.csv('GPL14951_c.csv', sep=',',  header=TRUE,
                        na.strings=c('','NA','---'))

GPL14951 <- rbind(GPL14951_a, GPL14951_b, GPL14951_c)



# merge the pancreatic files into one file for the series and platform
pancreatic <- merge(GSE131859, GPL14951, by.x='ID_REF', by.y='ID')
colnames(pancreatic)
# [1] "ID_REF"                           "MiaPACA2_EV"                     
# [3] "MiaPACA2_EV.Detection_Pval"       "MiaPACA2_ARHGEF10"               
# [5] "MiaPACA2_ARHGEF10.Detection_Pval" "Hs766T_shGFP"                    
# [7] "Hs766T_shGFP.Detection_Pval"      "Hs766T_shARHGEF10"               
# [9] "Hs766T_shARHGEF10.Detection_Pval" "X"                               
# [11] "Transcript"                       "Species"                         
# [13] "Source"                           "Search_Key"                      
# [15] "ILMN_Gene"                        "Source_Reference_ID"             
# [17] "RefSeq_ID"                        "Entrez_Gene_ID"                  
# [19] "GI"                               "Accession"                       
# [21] "Symbol"                           "Protein_Product"                 
# [23] "Array_Address_Id"                 "Probe_Type"                      
# [25] "Probe_Start"                      "SEQUENCE"                        
# [27] "Chromosome"                       "Probe_Chr_Orientation"           
# [29] "Probe_Coordinates"                "Cytoband"                        
# [31] "Definition"                       "Ontology_Component"              
# [33] "Ontology_Process"                 "Ontology_Function"               
# [35] "Synonyms"                         "Obsolete_Probe_Id"               
# [37] "GB_ACC"      

# keep only the Symbol, pancreatic cancer types for features
pancreatic1 <- pancreatic[,c(2:9,21)] # cell culture treated after extraction from pancreatic cncr
# the MiaPACA2_EV is empty vector control, MiaPACA2_ARHGEF10 is overexpressing, 
# Hs766T_shGFP is the short hairpin RNA (shRNA) knockdown of ARHGEF10 control (sep. cntrl smpl),
# and Hs766T_shARHGEF10 is the knockdown of ARHGEF10 

pancreatic1 <- pancreatic1[!duplicated(pancreatic1),]

library(dplyr)

pancr <- pancreatic1 %>% group_by(Symbol) %>% summarise(n=n())
pancr2 <- pancreatic1 %>% group_by(Symbol) %>% summarise_at(vars(MiaPACA2_EV:Hs766T_shARHGEF10.Detection_Pval),
                                                            mean)

pancr3 <- merge(pancr, pancr2, by.x='Symbol', by.y='Symbol')

write.csv(pancr3, 'pancreaticCancerGenes.csv', row.names=FALSE)

pancreatic4 <- pancr3 %>% mutate(DE_Incr_ARHGEF10_cntrl=MiaPACA2_ARHGEF10-MiaPACA2_EV)
pancreatic5 <- pancreatic4 %>% mutate(DE_decr_ARGHEF10_cntrl=Hs766T_shARHGEF10-Hs766T_shGFP)
pancreatic6 <- pancreatic5 %>% mutate(FC_incr_ARHGEF10_cntrl = MiaPACA2_ARHGEF10/MiaPACA2_EV)
pancreatic7 <- pancreatic6 %>% mutate(FC_decr_ARGHEF10_cntrl = Hs766T_shARHGEF10/Hs766T_shGFP)
row.names(pancreatic7) <- pancreatic7$Symbol
pancreatic7 <- pancreatic7[,c(11:14,3:10,1:2)]
write.csv(pancreatic7,'pancreaticCancerGeneStats.csv', row.names=TRUE)

#endoPancreatic <- merge(pancr3, EndocrineGenes, by.x='Symbol', by.y='Gene')
endopancrStats <- merge(pancreatic7, EndocrineGenes, by.x='Symbol', by.y='Gene')
row.names(endopancrStats) <- endopancrStats$Symbol
write.csv(endopancrStats, 'endocrinePancreaticCancerGeneStats.csv', row.names=TRUE)

#CBD_Pancreatic <- merge(painGenes, pancr3, by.x='Gene', by.y='Symbol')
CBD_pancreaticStats <- merge(pancreatic7, painGenes, by.x='Symbol', by.y='Gene')
row.names(CBD_pancreaticStats) <- CBD_pancreaticStats$Gene
write.csv(CBD_pancreaticStats, 'CBD_pancreatic_Gene_Stats.csv', row.names=TRUE)

# all genes in endocrine and CBD sets increased when ARHGEF10 increased, and decreased when 
# ARHGEF10 decreased with their respective control samples

# there are genes in the pancreatic Gene Stats set of all genes that incr when ARHGEF10 decr,
# and vice versa, and those that have no change.

library(ggplot2)
pancreaticOrder <- pancreatic7[order(pancreatic7$FC_incr_ARHGEF10_cntrl, decreasing=FALSE),]
pancreaticOrder1 <- subset(pancreaticOrder, pancreaticOrder$FC_incr_ARHGEF10_cntrl<0.3)

png('GenesNegCorrPancreaticCancerARHGEF10.png', width=768, height=576)

g <- ggplot(pancreaticOrder1, aes(x=FC_incr_ARHGEF10_cntrl, y= as.factor(Symbol)))
g= g+ ggtitle('Genes Negatively Correlated with Increased Tumor Suppressor ARHGEF10 
          in Cell Culture Pancreatic Cancer')
g= g+ ylab('Genes Negatively Correlated to Increased ARHGEF10')
g= g+ xlab('Fold Change Decreased More than 70 Percent')
g= g+ geom_point(aes(colour=DE_Incr_ARHGEF10_cntrl),size=6, alpha=0.9)
g

dev.off()

pancreaticOrder2 <- pancreatic7[order(pancreatic7$FC_incr_ARHGEF10_cntrl, decreasing=TRUE),]

pancreaticOrder3 <- subset(pancreaticOrder, pancreaticOrder$FC_incr_ARHGEF10_cntrl>50)


png('GenesPosCorrPancreaticCancerARHGEF10.png', width=768, height=576)

g <- ggplot(pancreaticOrder3, aes(x=FC_incr_ARHGEF10_cntrl, y= as.factor(Symbol)))
g= g+ ggtitle('Genes Positively Correlated with Increased Tumor Suppressor ARHGEF10 
          in Cell Culture Pancreatic Cancer')
g= g+ ylab('Genes Positively Correlated to Increased ARHGEF10')
g= g+ xlab('Fold Change Increased More than 50 Fold')
g= g+ geom_point(aes(colour=DE_Incr_ARHGEF10_cntrl),size=6, alpha=0.9)
g

dev.off()

CBD_pancreaticStats <- read.csv('CBD_pancreatic_Gene_Stats.csv', sep=',', header=TRUE,
                                row.names=1)
CBD_pancrStatsOrder <- CBD_pancreaticStats[order(CBD_pancreaticStats$FC_incr_ARHGEF10_cntrl,
                                                 decreasing=TRUE),]
CBD <- subset(CBD_pancrStatsOrder, CBD_pancrStatsOrder$FC_incr_ARHGEF10_cntrl > 5)

library(ggplot2)

png('CBD_GenesPosCorrPancreaticCancerARHGEF10.png', width=768, height=576)

g <- ggplot(CBD, aes(x=FC_incr_ARHGEF10_cntrl, y= as.factor(Symbol)))
g= g+ ggtitle('CBD Genes Positively Correlated with Increased Tumor Suppressor ARHGEF10 
          in Cell Culture Pancreatic Cancer')
g= g+ ylab('Genes Positively Correlated to Increased ARHGEF10')
g= g+ xlab('Fold Change Increased More than 5 Fold')
g= g+ geom_point(aes(colour=DE_Incr_ARHGEF10_cntrl),size=6, alpha=0.9)
g

dev.off()

endopancrStats <- read.csv('endocrinePancreaticCancerGeneStats.csv', sep=',', header=TRUE,
                           row.names=1)
endo <- endopancrStats[order(endopancrStats$FC_incr_ARHGEF10_cntrl, decreasing=TRUE),]
endo1 <- subset(endo, endo$FC_incr_ARHGEF10_cntrl>5)

png('Endocrine_GenesPosCorrPancreaticCancerARHGEF10.png', width=768, height=576)

g <- ggplot(endo1, aes(x=FC_incr_ARHGEF10_cntrl, y= as.factor(Symbol)))
g= g+ ggtitle('Endocrine Genes Positively Correlated with Increased Tumor Suppressor ARHGEF10 
          in Cell Culture Pancreatic Cancer')
g= g+ ylab('Genes Positively Correlated to Increased ARHGEF10')
g= g+ xlab('Fold Change Increased More than 5 Fold')
g= g+ geom_point(aes(colour=DE_Incr_ARHGEF10_cntrl),size=6, alpha=0.9)
g

dev.off()

###############################################################################################

# go back to control of pancreatic cancer genes
pancreatic7 <- read.csv('pancreaticCancerGeneStats.csv', row.names=1, sep=',', header=TRUE)

# pull in the CBD pain and endocrine gene datasets
painGenes <- read.csv('GenesSumm1.csv', header=TRUE, sep=',')
colnames(painGenes)[1] <- 'Gene'
EndocrineGenes <- read.csv('GenesSumm2.csv', header=TRUE, sep=',')
colnames(EndocrineGenes)[1] <- 'Gene'

# Also pull in the stats on the nonUL genes to get the mean of the normal nonUL female genes
# roughly compare the non-treated pancreatic cancer genes of one sample to the mean of 
# the non-uterine-leiomyoma genes of 51 females by mean of each gene. The pancreatic cancer genes
# are the row means of duplicate genes in original file, both are microarray samples

UL_stats <- read.csv('all_common_12173_130_fold_magnitude.csv', sep=',', 
                      row.names=1,header=TRUE)
row.names(UL_stats) <- UL_stats$GENE_SYMBOL
colnames(UL_stats)[1:10]
# [1] "GENE"        "CYTOBAND"    "GENE_SYMBOL" "Counts"      "UL_mean"     "nonUL_mean" 
# [7] "DE"          "Magnitude"   "foldChange"  "GSM1667144"

LeiomyomaStats <- UL_stats[,c(3,5:9)]
colnames(LeiomyomaStats)[c(4,6)] <- c('DE_UL','FoldChange_UL')
write.csv(LeiomyomaStats,'LeiomyomaStats.csv', row.names=TRUE)

# Also pull in the stats on glioma means from non glioma 
gliomaCBD <- read.csv('gliomaCBD.csv', row.names=1, header=TRUE, sep=',')

# Use the control means of glioma tumors not treated with CBD
gliomaControl <- gliomaCBD[,c(2,4,6)]
gliomaMeans <- data.frame(rowMeans(gliomaControl))
colnames(gliomaMeans) <- 'gliomaMeans'
gliomaMeans$Gene <- row.names(gliomaMeans)
gliomaMeans$Gene <- gsub(' ','', gliomaMeans$Gene)
gliomaMeans$Gene <- as.factor(gliomaMeans$Gene)

gliomaCBDtreated <- gliomaCBD[,c(1,3,5)]
gliomaCBDmeans <- data.frame(rowMeans(gliomaCBDtreated))
colnames(gliomaCBDmeans) <- 'gliomaMeans_CBD'
gliomaCBDmeans$Gene <- row.names(gliomaCBDmeans)
gliomaCBDmeans$Gene <- gsub(' ', '', gliomaCBDmeans$Gene)
gliomaCBDmeans$Gene <- as.factor(gliomaCBDmeans$Gene)

#use the control genes of pancreatic tumor values not incr/decr with ARHGEF10
pancreaticControl <- data.frame(pancreatic7[,5])
row.names(pancreaticControl) <- row.names(pancreatic7)
colnames(pancreaticControl) <- 'pancreaticMeans'
pancreaticControl$Gene <- row.names(pancreaticControl)
pancreaticControl$Gene <- as.factor(pancreaticControl$Gene)

# read in the sebum CBD data set and use the sebum controls not treated with CBD
sebumCBD <- read.csv('sebumCBD.csv', row.names=1, header=TRUE, sep=',')
sebumControl <- data.frame(sebumCBD[,4])
row.names(sebumControl) <- row.names(sebumCBD)
colnames(sebumControl) <- 'sebumMeans'
sebumControl$Gene <- row.names(sebumControl)
sebumControl$Gene <- as.factor(sebumControl$Gene)

sebumCBDmeans <- data.frame(sebumCBD[,5])
row.names(sebumCBDmeans) <- row.names(sebumCBD)
colnames(sebumCBDmeans) <- 'sebum_CBDtreated_Means'
sebumCBDmeans$Gene <- row.names(sebumCBDmeans)
sebumCBDmeans$Gene <- as.factor(sebumCBDmeans$Gene)

LeiomyomaControl <- data.frame(LeiomyomaStats[,c(1:3)])

# merge these datasets to compare non treated control samples from different microarray studies
controls <- merge(LeiomyomaControl, pancreaticControl, by.x='GENE_SYMBOL', by.y='Gene')
controls2 <- merge(controls, sebumControl, by.x='GENE_SYMBOL', by.y='Gene' )
controls3 <- merge(controls2, gliomaMeans, by.x='GENE_SYMBOL', by.y='Gene')

write.csv(controls3, 'controls_panc_sebm_UL_glma.csv', row.names=FALSE)

# Now we have a data set of different tissu and tumor control means
# add the other statistical data of CBD treated means for glioma and sebum
# gliomaCBDmeans and sebumCBDmeans datasets 
# and the values of adding ARHGEF10 --- known to suppress pancreatic tumors
# then compare the gene values, which incr/decr in samples and what tissues

controls4 <- merge(controls3, gliomaCBDmeans, by.x='GENE_SYMBOL', by.y='Gene')
controls5 <- merge(controls4, sebumCBDmeans, by.x='GENE_SYMBOL', by.y='Gene')

pancreaticStats <- pancreatic7[,c(1:4,13)]
colnames(pancreaticStats)[1:4] <- c('DE_incr_ARHGFE10_pancr','DE_decr_ARHGFE10_pancr',
                                    'FC_incr_ARHGFE10_pancr','FC_decr_ARHGFE10_pancr')

controls6 <- merge(controls5, pancreaticStats, by.x='GENE_SYMBOL', by.y='Symbol')

controls7 <- controls6[,c(1:3,6,7,5,8,4,9:12)]
colnames(controls7)
# [1] "GENE_SYMBOL"            "UL_mean"                "nonUL_mean"            
# [4] "gliomaMeans"            "gliomaMeans_CBD"        "sebumMeans"            
# [7] "sebum_CBDtreated_Means" "pancreaticMeans"        "DE_incr_ARHGFE10_pancr"
# [10] "DE_decr_ARHGFE10_pancr" "FC_incr_ARHGFE10_pancr" "FC_decr_ARHGFE10_pancr"

write.csv(controls7, 'UL_glioma_sebm_pancr_stats.csv', row.names=FALSE)#9281X12


#######################################################################################

# look at the breast cancer, colon cancer, and stomach cancer array data
# import them and subset each into datasets less than 25 MB for github
# then import in with each subset to run code

# breast cancer array data for the series and platform from the gene expression omnibus, GEO
# GPL571-17391.txt and GSE110332_series_matrix.txt search GPL571 and GSE110332 in GEO

GSE110332 <- read.delim('GSE110332_series_matrix.txt',sep='\t', quote="", header=TRUE,
                        comment.char='!', na.strings=c('','NA'))
colnames(GSE110332) <- gsub('X.','', colnames(GSE110332))
colnames(GSE110332) <- gsub('.','', colnames(GSE110332), fixed=TRUE)
GSE110332$ID_REF <- gsub('"','',as.character(GSE110332$ID_REF), fixed=TRUE, perl=TRUE)

GPL571 <- read.delim('GPL571-17391.txt', sep='\t', quote="", header=TRUE,
                       comment.char='#', na.strings=c('','NA','---'))

GPL571_a <- GPL571[1:7999,]
GPL571_a_2 <- GPL571[8000:11139,]
GPL571_b <- GPL571[11140:22277,]

write.csv(GPL571_a, 'GPL571_a.csv', row.names=FALSE)
write.csv(GPL571_a_2, 'GPL571_a_2.csv', row.names=FALSE)
write.csv(GPL571_b, 'GPL571_b.csv', row.names=FALSE)

GPL571_a <- read.csv('GPL571_a.csv', sep=',', header=TRUE)
GPL571_a_2 <- read.csv('GPL571_a_2.csv', sep=',', header=TRUE)
GPL571_b <- read.csv('GPL571_b.csv', sep=',', header=TRUE)

brstCncr <- rbind(GPL571_a, GPL_a_2, GPL571_b)
brstCancer <- merge(GSE110332, brstCncr, by.x='ID_REF', by.y='ID')

write.csv(brstCancer, 'breastCancer.csv', row.names=FALSE)

# Stomach cancer array data, GPL13497-9755.txt and GSE64916_series_matrix.txt
GPL13497 <- read.delim('GPL13497-9755.txt', sep='\t', quote="", header=TRUE,
                       comment.char='#', na.strings=c('','NA','---'))
GSE64916 <- read.delim('GSE64916_series_matrix.txt', sep='\t', quote='', header=TRUE,
                       comment.char='!', na.strings=c('','NA','---'))
colnames(GSE64916) <- gsub('X.','', colnames(GSE64916))
colnames(GSE64916) <- gsub('.','', colnames(GSE64916), fixed=TRUE)
GSE64916$ID_REF <- gsub('"','', as.character(GSE64916$ID_REF), fixed=TRUE, perl=TRUE) 

stomachCancer <- merge(GSE64916, GPL13497, by.x='ID_REF', by.y='ID')

write.csv(stomachCancer,'StomachCancer.csv', row.names=FALSE)  

# colon cancer array data
GSE135749 <- read.delim('GSE135749_series_matrix.txt', sep='\t', quote='', header=TRUE,
                       comment.char='!', na.strings=c('','NA','---'))
colnames(GSE135749) <- gsub('X.','', colnames(GSE135749))
colnames(GSE135749) <- gsub('.','', colnames(GSE135749), fixed=TRUE)
GSE135749$ID_REF <- gsub('"','', as.character(GSE135749$ID_REF), fixed=TRUE, perl=TRUE) 

GPL10558 <- read.delim('GPL10558-50081.txt', sep='\t', quote="", header=TRUE,
                       comment.char='#', na.strings=c('','NA','---'))
#divide GPL10558 into 3 parts a:c bc it is 68.8 MB

GPL10558_a <- GPL10558[1:10000,]
GPL10558_a_2 <- GPL10558[10001:16036,]

GPL10558_b <- GPL10558[16037:25000,]
GPL10558_b_2 <- GPL10558[25001:32073,]

GPL10558_c <- GPL10558[32074:48107,]

write.csv(GPL10558_a, 'GPL10558_a.csv', row.names=FALSE)
write.csv(GPL10558_a_2, 'GPL10558_a_2.csv', row.names=FALSE)

write.csv(GPL10558_b, 'GPL10558_b.csv', row.names=FALSE)
write.csv(GPL10558_b_2, 'GPL10558_b_2.csv', row.names=FALSE)

write.csv(GPL10558_c, 'GPL10558_c.csv', row.names=FALSE)

GPL10558_a <- read.csv('GPL10558_a.csv', sep=',', header=TRUE)
GPL10558_a_2 <- read.csv('GPL10558_a_2.csv', sep=',', header=TRUE)

GPL10558_b <- read.csv('GPL10558_b.csv', sep=',', header=TRUE)
GPL10558_b_2 <- read.csv('GPL10558_b_2.csv', sep=',', header=TRUE)

GPL10558_c <- read.csv('GPL10558_c.csv', sep=',', header=TRUE)

GPL10558 <- rbind(GPL10558_a, GPL10558_a_2, GPL10558_b, GPL10558_b_2, GPL10558_c)

colonCancer <- merge(GSE135749, GPL10558, by.x='ID_REF', by.y='ID')

write.csv(colonCancer, 'ColonCancer.csv', row.names=FALSE)

#######################################################################################

# 'UL_glioma_sebm_pancr_stats.csv' and the newer files 'ColonCancer.csv' , 'breastCancer.csv' ,
# and 'StomachCancer.csv' are to be cleaned and combined to compare expression values of genes
# and their changes relative to the tissue type sampled and the tumor sampled in array data

UL_glm_sbm_pncr_stats <- read.csv('UL_glioma_sebm_pancr_stats.csv', sep=',', header=TRUE)

breastCancer <- read.csv('breastCancer.csv', sep=',', header=TRUE)
colonCancer <- read.csv('ColonCancer.csv', sep=',', header=TRUE)
stomachCancer <- read.csv('StomachCancer.csv', sep=',', header=TRUE)

brcr <- breastCancer[na.omit(breastCancer$Gene.Symbol),]
brscr <- brcr[,c(2:7,17)]
brstCncr <- brscr[!duplicated(brscr$Gene.Symbol),]
Gene <- strsplit(as.character(brstCncr$Gene.Symbol), '///', perl=TRUE)
gene <- lapply(Gene, '[',1)
brstCncr$Gene <- t(data.frame(gene))
breastCancer <- brstCncr[,-7]

write.csv(breastCancer,'BreastCancerCleaned.csv', row.names=FALSE)

ColonCancer <- colonCancer[,c(2:7,19)]

write.csv(ColonCancer,'ColonCancerCleaned.csv', row.names=FALSE)

Stomach <- stomachCancer[na.omit(stomachCancer$GENE_SYMBOL),]
StomachCancer <- Stomach[,c(2:6,12)]
stomachCncr <- StomachCancer[complete.cases(StomachCancer),]

write.csv(stomachCncr,'StomachCancerCleaned.csv', row.names=FALSE)


# the breast cancer data set 'BreastCancerCleaned.csv' has SH3GL2 as a tumor suppressor
# in the last 3 samples, and controls as the first 3 samples

breastCancer <- read.csv('BreastCancerCleaned.csv',sep=',', header=TRUE, na.strings=c('','NA'))
breastCancer <- breastCancer[complete.cases(breastCancer),]
brcr_cntrl <- breastCancer[,c(1:3,7)]
brcr_SH3GL2 <- breastCancer[,c(4:7)]

row.names(brcr_cntrl) <- brcr_cntrl$Gene
library(dplyr)

br_cr <- brcr_cntrl %>% group_by(Gene) %>% summarise_at(vars(GSM2987594:GSM2987596),
                                                        mean, na.rm=TRUE)
br_ctrl <- data.frame(br_cr)
row.names(br_ctrl) <- br_cr$Gene
br_ctrl <- br_ctrl[,-1]
br_ctrl$BRCNCR_CTRL_Mean <- rowMeans(br_ctrl)
br_ctrl$Gene <- row.names(br_ctrl)
br_ctrl <- br_ctrl[,-c(1:3)]


br_sh3 <- brcr_SH3GL2 %>% group_by(Gene) %>% summarise_at(vars(GSM2987597:GSM2987599),
                                                          mean, na.rm=TRUE)
br_cr_SH3GL2 <- data.frame(br_sh3)
row.names(br_cr_SH3GL2) <- br_cr_SH3GL2$Gene
br_cr_SH3GL2 <- br_cr_SH3GL2[,-1]
br_cr_SH3GL2$BRCNCR_SH3GL2_Mean <- rowMeans(br_cr_SH3GL2)
br_cr_SH3GL2$Gene <- row.names(br_cr_SH3GL2)
br_cr_SH3GL2 <- br_cr_SH3GL2[,-c(1:3)]

BRCR <- merge(br_ctrl, br_cr_SH3GL2, by.x='Gene', by.y='Gene')

UL_glm_sbm_pncr_stats <- read.csv('UL_glioma_sebm_pancr_stats.csv', sep=',', header=TRUE)

Means_diseases_Genes <- merge(BRCR, UL_glm_sbm_pncr_stats, by.x='Gene', by.y='GENE_SYMBOL')

# clean the colon cancer data
colonCancer <- read.csv('ColonCancer.csv', sep=',', header=TRUE)
# first two samples are controls, next two are short hair pin knockdowns of LGR5-1,
# last two are short hair pin RNA knockdowns of LGR5-2, this study gene changes, hence the
# negative expression values, that are also quantile normalized from
# this knockdown of two types to the LoVo colon cancer cell line

colonCancer <- colonCancer[,c(2:7,19)]
colonCancer <- colonCancer[!duplicated(colonCancer),]
colonCancer <- colonCancer[complete.cases(colonCancer),]
cln_cncr <- colonCancer %>% group_by(Symbol) %>% summarise_at(vars(GSM4027437:GSM4027446),
                                                              mean, na.rm=TRUE)
cln_cncr <- data.frame(cln_cncr)
row.names(cln_cncr) <- cln_cncr$Symbol
cln_cncr <- cln_cncr[,-1]

cln_cncr_ctrl <- cln_cncr[,1:2]
cln_cncr_ctrl$colon_cncr_ctrl_Mean <- rowMeans(cln_cncr_ctrl)
cln_cncr_ctrl$Gene <- row.names(cln_cncr_ctrl)
cln_cncr_ctrl <- cln_cncr_ctrl[,-c(1:2)]

cln_cncr_LGR51_loss <- cln_cncr[,3:4]
cln_cncr_LGR51_loss$colon_cncr_LGR51_loss_Mean <- rowMeans(cln_cncr_LGR51_loss)
cln_cncr_LGR51_loss$Gene <- row.names(cln_cncr_LGR51_loss)
cln_cncr_LGR51_loss <- cln_cncr_LGR51_loss[,-c(1:2)]

cln_cncr_LGR52_loss <- cln_cncr[,5:6]
cln_cncr_LGR52_loss$colon_cncr_LGR52_loss_Mean <- rowMeans(cln_cncr_LGR52_loss)
cln_cncr_LGR52_loss$Gene <- row.names(cln_cncr_LGR52_loss)
cln_cncr_LGR52_loss <- cln_cncr_LGR52_loss[,-c(1:2)]

cln <- merge(cln_cncr_ctrl, cln_cncr_LGR51_loss, by.x='Gene', by.y='Gene')
cln2 <- merge(cln, cln_cncr_LGR52_loss, by.x='Gene', by.y='Gene')

diseaseGenes <- merge(Means_diseases_Genes,cln2, by.x='Gene', by.y='Gene')


# Now add the stomach cancer array data
# the first 4 samples are gastric stomach cancer genes in the peripheral blood, 
# the last, sample 5, is a healthy control, whole blood samples, ages 39-46 for all 5
stomachCancer <- read.csv('StomachCancerCleaned.csv', sep=',', 
                          header=TRUE, na.strings=c('','NA'))
Stomach <- stomachCancer[complete.cases(stomachCancer),]
Stomach_ctrl <- Stomach[,5:6]
colnames(Stomach_ctrl) <- c('Stomach_Blood_Ctrl_Mean','Gene')

Stomach_cncr <- Stomach[,c(1:4,6)]
stmch_cncr <- Stomach_cncr %>% group_by(GENE_SYMBOL) %>% summarise_at(vars(GSM1583284:GSM1583287), 
                                                                      mean, na.rm=TRUE)
Stomach_cncr <- data.frame(stmch_cncr)
row.names(Stomach_cncr) <- Stomach_cncr$GENE_SYMBOL
Stomach_cncr <- Stomach_cncr[,-1]
Stomach_cncr$Stomach_Cancer_Mean <- rowMeans(Stomach_cncr)
Stomach_cncr$Gene <- row.names(Stomach_cncr)
Stomach_cncr <- Stomach_cncr[,-c(1:4)]
colnames(Stomach_cncr)[1] <- 'Stomach_Cancer_Blood_Mean'

Stomach <- merge(Stomach_cncr, Stomach_ctrl, by.x='Gene', by.y='Gene')

DiseaseGenes <- merge(diseaseGenes, Stomach, by.x='Gene', by.y='Gene')

write.csv(DiseaseGenes, 'br_ul_panc_stmc_cln_glm_sbm_Means.csv', row.names=FALSE)
# ARHGEF10 is a tumor suppressor in pancreatic cancer tumors, with incr and decr samples
# most genes also increase and decrease with the incr or decr in ARHGEF10, some move opposite
# LGR5 is shown in colon cancer for changes, The colon cancer data is of blood gene changes 
# using lentivirus infection in cell cultures from the colon cancer blood, 
# breast cancer samples are compared to SH3GL2 as a tumor supporessor in BRCR
# the glioma and sebum samples are compared to CBD treated samples of same
# glioma from the tumor biopsies and sebum from hair follicle extractions of sebacious glands
# for acne experiment using CBD treated sebum
# the UL are uterine leiomyomas, and the nonUL are otherwise healthy samples, both from uterine
# tissue samples

###################################################################################################
###################################################################################################
###################################################################################################


# read in all samples from all data sets and omit the statistical data added like fold change
# differential expression, means, and magnitude, as well as counts when taking row means
# of duplicate genes in a data set

gliomaCBD <- read.csv('gliomaCBD.csv', row.names=1, header=TRUE, sep=',')
sebumCBD <- read.csv('sebumCBD.csv', row.names=1, header=TRUE, sep=',')
UL_FC <- read.csv('all_common_12173_130_fold_magnitude.csv', sep=',', header=TRUE, row.names=1)
pancr3 <- read.csv('pancreaticCancerGenes.csv', sep=',', header=TRUE, na.strings=c('','NA'))
brstCancer <- read.csv('BreastCancerCleaned.csv', sep=',', header=TRUE, na.strings=c('', 'NA'))
clnCancer <- read.csv('ColonCancerCleaned.csv', sep=',', header=TRUE, na.strings=c('','NA'))
stmCancer <- read.csv('StomachCancerCleaned.csv', sep=',', header=TRUE, na.strings=c('','NA'))

gliomaCBD$gene <-row.names(gliomaCBD)
gliomaCBD$gene <- gsub(' ','', as.character(gliomaCBD$gene))
gliomaCBD$gene <- as.factor(gliomaCBD$gene)

sebumCBD$gene <- row.names(sebumCBD)
sebumCBD$gene <- as.factor(sebumCBD$gene)

colnames(gliomaCBD) <- c("GSM1399027.glioma.CBD" , "GSM1399028.glioma.CTRL", 
                         "GSM1399029.glioma.CBD" , "GSM1399030.glioma.CTRL" ,
                         "GSM1399031.glioma.CBD","GSM1399032.glioma.CTRL", "gene" )

colnames(sebumCBD) <- c("GSM1385013.sebum.REP1", "GSM1385014.sebum.REP2", "GSM1385015.sebum.REP2",
                        "GSM1385016.sebum.CTRL", "GSM1385017.sebum.CBD", "gene")

colnames(clnCancer) <- c("GSM4027437.colon.ctrl", "GSM4027439.colon.ctrl", 
                         "GSM4027441.colon.LGR5.2.down", "GSM4027442.colon.LGR5.2.down", 
                         "GSM4027444.colon.LGR5.1.down", "GSM4027446.colon.LGR5.1.down", "Symbol")

colnames(pancr3) <- c("Symbol"  ,                         "n"      ,                         
                      "pancreas_ctrl"  ,                    "pancreas_ctrl_Pval" ,     
                      "pancreas_up_ARHGEF10",                "pancreas_up_ARHGEF10_Pval",
                      "pancreas_up_then_down_ARHGEF10"  , "pancreas_up_then_down_ARHGEF10_Pval" ,    
                      "pancreas_down_ARHGEF10" ,               "pancreas_down_ARHGEF10_Pval")

#remove pvalues from pancreatic data frame
pancr3 <- pancr3[,c(1,3,5,7,9)]

colnames(stmCancer) <- c("GSM1583284.stomach.cancer",  "GSM1583285.stomach.cancer" , 
                         "GSM1583286.stomach.cancer" , "GSM1583287.stomach.cancer",  
                         "GSM1583288.stomach.healthy" , "GENE_SYMBOL")

colnames(brstCancer) <- c("GSM2987594.ctrl.brstCncr", "GSM2987595.ctrl.brstCncr", 
                          "GSM2987596.ctrl.brstCncr" ,
                          "GSM2987597.SH3GL2.treated.brstCncr", 
                          "GSM2987598.SH3GL2.treated.brstCncr" ,
                          "GSM2987599.SH3GL2.treated.brstCncr" ,"Gene" )

# remove meta fields from UL dataset of stats
UL <- UL_FC[,c(3,10:130)]

colnames(UL)[2:52] <- paste(as.character(colnames(UL)[2:52]), sep='_', 'nonUL')
colnames(UL)[53:122] <- gsub('UL','',as.character(colnames(UL)[53:122]))
colnames(UL)[53:122] <- paste(as.character(colnames(UL)[53:122]), sep='_', 'UL')

# Now that all the samples are identified, the colon cancer and the sebum values are change values,
# hence the negative values, they will be added as they are

# merge all the datasets into one big set by their respective gene symbols

df <- merge(UL, stmCancer, by.x='GENE_SYMBOL', by.y='GENE_SYMBOL')
df1 <- merge(df, brstCancer, by.x='GENE_SYMBOL', by.y='Gene')
df2 <- merge(df1, pancr3, by.x='GENE_SYMBOL', by.y='Symbol')
df3 <- merge(df2, gliomaCBD, by.x='GENE_SYMBOL', by.y='gene')
df4 <- merge(df3, clnCancer, by.x='GENE_SYMBOL', by.y='Symbol')
df5 <- merge(df4, sebumCBD, by.x='GENE_SYMBOL', by.y='gene')

# there were duplicate entries per gene in each sample to add the additional samples,

#remove duplicate gene entries
df5 <- df5[!duplicated(df5$GENE_SYMBOL),]
row.names(df5) <- df5$GENE_SYMBOL

#remove the GENE_SYMBOL field as the rows now identify the gene of each sample
all <- df5[,-c(1)]

#write this dataset out to csv
write.csv(all, 'all_samples_5871_153.csv', row.names=TRUE)

###################################################################################################
###################################################################################################
###################################################################################################

all_samples <- read.csv('all_samples_5871_153.csv', row.names=1)

colnames(all_samples)
# [1] "GSM1667144_nonUL"                   "GSM1667145_nonUL"                  
# [3] "GSM1667146_nonUL"                   "GSM336252_nonUL"                   
# [5] "GSM336253_nonUL"                    "GSM336254_nonUL"                   
# [7] "GSM336255_nonUL"                    "GSM336256_nonUL"                   
# [9] "GSM336257_nonUL"                    "GSM336258_nonUL"                   
# [11] "GSM336259_nonUL"                    "GSM336260_nonUL"                   
# [13] "GSM336261_nonUL"                    "GSM336262_nonUL"                   
# [15] "GSM336263_nonUL"                    "GSM336264_nonUL"                   
# [17] "GSM336265_nonUL"                    "GSM336266_nonUL"                   
# [19] "GSM336267_nonUL"                    "GSM336268_nonUL"                   
# [21] "GSM336269_nonUL"                    "GSM336270_nonUL"                   
# [23] "GSM336271_nonUL"                    "GSM336272_nonUL"                   
# [25] "GSM336273_nonUL"                    "GSM336274_nonUL"                   
# [27] "GSM336275_nonUL"                    "GSM336276_nonUL"                   
# [29] "GSM336277_nonUL"                    "GSM336278_nonUL"                   
# [31] "GSM52661_nonUL"                     "GSM52662_nonUL"                    
# [33] "GSM52663_nonUL"                     "GSM52664_nonUL"                    
# [35] "GSM52665_nonUL"                     "GSM52666_nonUL"                    
# [37] "GSM52667_nonUL"                     "GSM52668_nonUL"                    
# [39] "GSM52669_nonUL"                     "GSM52670_nonUL"                    
# [41] "GSM52671_nonUL"                     "GSM9098_nonUL"                     
# [43] "GSM9099_nonUL"                      "GSM9100_nonUL"                     
# [45] "GSM9101_nonUL"                      "GSM9102_nonUL"                     
# [47] "GSM569424_nonUL"                    "GSM569425_nonUL"                   
# [49] "GSM569426_nonUL"                    "GSM569427_nonUL"                   
# [51] "GSM569428_nonUL"                    "GSM1667147_UL"                     
# [53] "GSM1667148_UL"                      "GSM1667149_UL"                     
# [55] "GSM336202_UL"                       "GSM336203_UL"                      
# [57] "GSM336204_UL"                       "GSM336205_UL"                      
# [59] "GSM336206_UL"                       "GSM336207_UL"                      
# [61] "GSM336208_UL"                       "GSM336209_UL"                      
# [63] "GSM336210_UL"                       "GSM336211_UL"                      
# [65] "GSM336212_UL"                       "GSM336213_UL"                      
# [67] "GSM336214_UL"                       "GSM336215_UL"                      
# [69] "GSM336216_UL"                       "GSM336217_UL"                      
# [71] "GSM336218_UL"                       "GSM336219_UL"                      
# [73] "GSM336220_UL"                       "GSM336221_UL"                      
# [75] "GSM336222_UL"                       "GSM336223_UL"                      
# [77] "GSM336224_UL"                       "GSM336225_UL"                      
# [79] "GSM336226_UL"                       "GSM336227_UL"                      
# [81] "GSM336228_UL"                       "GSM336229_UL"                      
# [83] "GSM336230_UL"                       "GSM336231_UL"                      
# [85] "GSM336232_UL"                       "GSM336233_UL"                      
# [87] "GSM336234_UL"                       "GSM336235_UL"                      
# [89] "GSM336236_UL"                       "GSM336237_UL"                      
# [91] "GSM336238_UL"                       "GSM336239_UL"                      
# [93] "GSM336240_UL"                       "GSM336241_UL"                      
# [95] "GSM336242_UL"                       "GSM336243_UL"                      
# [97] "GSM336244_UL"                       "GSM336245_UL"                      
# [99] "GSM336246_UL"                       "GSM336247_UL"                      
# [101] "GSM336248_UL"                       "GSM336249_UL"                      
# [103] "GSM336250_UL"                       "GSM336251_UL"                      
# [105] "GSM38689_UL"                        "GSM38690_UL"                       
# [107] "GSM38691_UL"                        "GSM38692_UL"                       
# [109] "GSM38693_UL"                        "GSM38694_UL"                       
# [111] "GSM38695_UL"                        "GSM9093_UL"                        
# [113] "GSM9094_UL"                         "GSM9095_UL"                        
# [115] "GSM9096_UL"                         "GSM9097_UL"                        
# [117] "GSM569429_UL"                       "GSM569430_UL"                      
# [119] "GSM569431_UL"                       "GSM569432_UL"                      
# [121] "GSM569433_UL"                       "GSM1583284.stomach.cancer"         
# [123] "GSM1583285.stomach.cancer"          "GSM1583286.stomach.cancer"         
# [125] "GSM1583287.stomach.cancer"          "GSM1583288.stomach.healthy"        
# [127] "GSM2987594.ctrl.brstCncr"           "GSM2987595.ctrl.brstCncr"          
# [129] "GSM2987596.ctrl.brstCncr"           "GSM2987597.SH3GL2.treated.brstCncr"
# [131] "GSM2987598.SH3GL2.treated.brstCncr" "GSM2987599.SH3GL2.treated.brstCncr"
# [133] "pancreas_ctrl"                      "pancreas_up_ARHGEF10"              
# [135] "pancreas_up_then_down_ARHGEF10"     "pancreas_down_ARHGEF10"            
# [137] "GSM1399027.glioma.CBD"              "GSM1399028.glioma.CTRL"            
# [139] "GSM1399029.glioma.CBD"              "GSM1399030.glioma.CTRL"            
# [141] "GSM1399031.glioma.CBD"              "GSM1399032.glioma.CTRL"            
# [143] "GSM4027437.colon.ctrl"              "GSM4027439.colon.ctrl"             
# [145] "GSM4027441.colon.LGR5.2.down"       "GSM4027442.colon.LGR5.2.down"      
# [147] "GSM4027444.colon.LGR5.1.down"       "GSM4027446.colon.LGR5.1.down"      
# [149] "GSM1385013.sebum.REP1"              "GSM1385014.sebum.REP2"             
# [151] "GSM1385015.sebum.REP2"              "GSM1385016.sebum.CTRL"             
# [153] "GSM1385017.sebum.CBD" 
