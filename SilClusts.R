library(tidyverse)
library(RColorBrewer)
library(cluster)

AllZT <- read.table("lab/Heatmaps_CV2_expression/BrenneckeMod_allZT_RUV2.txt",
                    header=TRUE, sep='\t')

AllZT_BioVar <- select(AllZT, Gene, contains("BioVar"))

AllZT_BioVar[is.na(AllZT_BioVar)]<-0

HVG <- read.table("lab/Heatmaps_CV2_expression/HVG_allZT_nbZT.txt",header=TRUE,sep='\t')

HVG_allZT_BioVar <- inner_join(AllZT_BioVar,HVG, by=c("Gene"="allZT")) %>%
  select(., Gene, contains("BioVar"))

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
#---------------------------------k-means----------------------------------------
x <- HVG_allZT_BioVar[2:13]

k_means <- kmeans(x, 4, iter.max = 10, nstart = 1)

#-------------------------------Silhouette--------------------------------------

sil <- silhouette(k_means$cluster,dist(HVG_allZT_BioVar[2:13]))
windows()
plot(sil)

silhouette_score <- function(k){
  km <- kmeans(HVG_allZT_BioVar[2:13], centers = k, nstart=25)
  ss <- silhouette(km$cluster, dist(HVG_allZT_BioVar[2:13]))
  mean(ss[, 3])
}
k <- 2:10
avg_sil <- sapply(k, silhouette_score)
plot(k, type='b', avg_sil, xlab='Number of clusters', ylab='Average Silhouette Scores', frame=FALSE)

