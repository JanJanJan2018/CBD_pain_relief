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
