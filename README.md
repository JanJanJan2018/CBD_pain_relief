# CBD_pain_relief
research on CBD topical and oral treatment of pain and other well-being diseases

Gene Expression Omnibus data used and articles available from open source research sources also used in this repository to study the current literature on CBD therapies

GPL6244 is a GEO platform that was too large 94mb to upload to this repository. All GEO datasets are retrieved with GSE or GPL ID, download the .gz extensions of the txt files for the series and the platforms to access. This research used WinZip to unzip the gz files

The larger data set GSE56978 was cleaned of NAs and duplicates then split into four sections of 5538 rows each as csv files with the modified field extensions to what type of sample it is. The row names are false. Read in each of gse56978_1.csv, gse56978_2.csv, gse56978_3.csv, and gse56978_4.csv then combine with '>GSE56978 <- rbind(gse56978_1, gse56978_2, gse56978_3,gse56978_4) assuming they are stored as data tables with those respective names. The entire file is >90MB, and this version of github only allows <=25mb in file size.

