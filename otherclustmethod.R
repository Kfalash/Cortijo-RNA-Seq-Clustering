#Other clustering look.

library(tidyverse)
library(ggplot2)
library(fpc)

#Original data
AllZT <- read.table("lab/Heatmaps_CV2_expression/BrenneckeMod_allZT_RUV2.txt",
                    header=TRUE, sep='\t')

AllZT_BioVar <- select(AllZT, Gene, contains("BioVar"))

AllZT_BioVar[is.na(AllZT_BioVar)]<-0

HVG <- read.table("lab/Heatmaps_CV2_expression/HVG_allZT_nbZT.txt",header=TRUE,sep='\t')

HVG_allZT_BioVar <- inner_join(AllZT_BioVar,HVG, by=c("Gene"="allZT")) %>%
  select(., Gene, contains("BioVar"))

#----------------------------Analysis------------------------------------------

#scale the variables
scaled_wd <- scale(HVG_allZT_BioVar[2:13])

#Hierarchical Clustering
d <- dist(scaled_wd, method = "euclidean")

#distance matrix
h_clust <- hclust(d, method = "average") #clustering

#dendrogram
plot(h_clust,labels = HVG_allZT_BioVar$Gene)

rect.hclust(h_clust,k=2)

#extract clusters
groups <- cutree(h_clust,k=2)

#pca
pcmp <- princomp(scaled_wd)
pred_pc <- predict(pcmp, newdata=scaled_wd)[,1:2]

comp_dt <- cbind(as.data.table(pred_pc),cluster = as.factor(groups), Labels = HVG_allZT_BioVar$Gene)
ggplot(comp_dt,aes(Comp.1,Comp.2))+ geom_point(aes(color = cluster),size=3)

#kmeans
kclust <- kmeans(scaled_wd,centers = 2,iter.max = 100)
ggplot(comp_dt,aes(Comp.1,Comp.2))+ geom_point(aes(color = as.factor(kclust$cluster)),size=3)

tunek <- kmeansruns(scaled_wd,krange = 1:10,criterion = "ch")
tunek$bestk #2
tunekw <- kmeansruns(scaled_wd,krange = 1:10,criterion = "asw")
tunekw$bestk #2
