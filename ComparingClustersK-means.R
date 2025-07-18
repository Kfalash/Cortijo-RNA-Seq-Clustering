library(tidyverse)
library(RColorBrewer)

AllZT <- read.table("lab/Heatmaps_CV2_expression/BrenneckeMod_allZT_RUV2.txt",
                    header=TRUE, sep='\t')

AllZT_BioVar <- select(AllZT, Gene, contains("BioVar"))

AllZT_BioVar[is.na(AllZT_BioVar)]<-0

HVG <- read.table("lab/Heatmaps_CV2_expression/HVG_allZT_nbZT.txt",header=TRUE,sep='\t')

HVG_allZT_BioVar <- inner_join(AllZT_BioVar,HVG, by=c("Gene"="allZT")) %>%
  select(., Gene, contains("BioVar"))

#---------------------------------k-means----------------------------------------
x <- HVG_allZT_BioVar[2:13]

k_means <- kmeans(x, 4, iter.max = 10, nstart = 1)

# ---------------------------Complete Clustering---------------------------------

HVG_AllZT_BioVar.pearson <- cor(t(HVG_allZT_BioVar[,c(2:13)]), method = "pearson")
HVG_AllZT_BioVar.pearson[is.na(HVG_AllZT_BioVar.pearson)]<-0
HVG_AllZT_BioVar.pearson.dist<-as.dist((1-HVG_AllZT_BioVar.pearson))
HVG_AllZT_BioVar.pearson.dist[is.na(HVG_AllZT_BioVar.pearson.dist)]<-0
HVG_AllZT_BioVar.pearson.dist.clust.complete <- hclust(HVG_AllZT_BioVar.pearson.dist, 
                                                       method="complete")
HVG_AllZT_BioVar.pearson.dist.clust.complete.den <- as.dendrogram(HVG_AllZT_BioVar.pearson.dist.clust.complete)
HVG_AllZT_BioVar.pearson.dist.clust.complete.col=brewer.pal(12, 'Set3')[cutree(HVG_AllZT_BioVar.pearson.dist.clust.complete, k = 4)]
HVG_allZT_BioVar$cluster_4_pearson <- cutree(HVG_AllZT_BioVar.pearson.dist.clust.complete, 
                                             k = 4)

Complete_Clusts <- HVG_allZT_BioVar[14]

# ----------------------------Average Clustering--------------------------------

HVG_AllZT_BioVar.pearson <- cor(t(HVG_allZT_BioVar[,c(2:13)]), method = "pearson")
HVG_AllZT_BioVar.pearson[is.na(HVG_AllZT_BioVar.pearson)]<-0
HVG_AllZT_BioVar.pearson.dist<-as.dist((1-HVG_AllZT_BioVar.pearson))
HVG_AllZT_BioVar.pearson.dist[is.na(HVG_AllZT_BioVar.pearson.dist)]<-0
HVG_AllZT_BioVar.pearson.dist.clust.complete <- hclust(HVG_AllZT_BioVar.pearson.dist, 
                                                       method="average")
HVG_AllZT_BioVar.pearson.dist.clust.complete.den <- as.dendrogram(HVG_AllZT_BioVar.pearson.dist.clust.complete)
HVG_AllZT_BioVar.pearson.dist.clust.complete.col=brewer.pal(12, 'Set3')[cutree(HVG_AllZT_BioVar.pearson.dist.clust.complete, k = 4)]
HVG_allZT_BioVar$cluster_4_pearson <- cutree(HVG_AllZT_BioVar.pearson.dist.clust.complete, 
                                             k = 4)

Average_Clusts <- HVG_allZT_BioVar[14]

# ------------------------------Ward Clustering---------------------------------

HVG_AllZT_BioVar.pearson <- cor(t(HVG_allZT_BioVar[,c(2:13)]), method = "pearson")
HVG_AllZT_BioVar.pearson[is.na(HVG_AllZT_BioVar.pearson)]<-0
HVG_AllZT_BioVar.pearson.dist<-as.dist((1-HVG_AllZT_BioVar.pearson))
HVG_AllZT_BioVar.pearson.dist[is.na(HVG_AllZT_BioVar.pearson.dist)]<-0
HVG_AllZT_BioVar.pearson.dist.clust.complete <- hclust(HVG_AllZT_BioVar.pearson.dist, 
                                                       method="ward.D")
HVG_AllZT_BioVar.pearson.dist.clust.complete.den <- as.dendrogram(HVG_AllZT_BioVar.pearson.dist.clust.complete)
HVG_AllZT_BioVar.pearson.dist.clust.complete.col=brewer.pal(12, 'Set3')[cutree(HVG_AllZT_BioVar.pearson.dist.clust.complete, k = 4)]
HVG_allZT_BioVar$cluster_4_pearson <- cutree(HVG_AllZT_BioVar.pearson.dist.clust.complete, 
                                             k = 4)

Ward_Clusts <- HVG_allZT_BioVar[14]

#--------------Intercetion of clusters with k-means and other clustering--------
#--------k-means and Average-------
table(k_means$cluster,Average_Clusts$cluster)

#-------k-means and Ward------
table(k_means$cluster,Ward_Clusts$cluster)

#-------k-means and Complete------
table(k_means$cluster,Complete_Clusts$cluster)
