
library(tidyverse)
library(factoextra)
library(clustertend)

### Clustering and heatmap of Corrected CV2 for HVGs ----

# Data import and formatting

df <- read.table("lab/Heatmaps_CV2_expression/BrenneckeMod_allZT_RUV2.txt",
                    header=TRUE, sep='\t')

AllZT_BioVar <- select(AllZT, Gene, contains("BioVar"))

AllZT_BioVar[is.na(AllZT_BioVar)]<-0

HVG <- read.table("lab/Heatmaps_CV2_expression/HVG_allZT_nbZT.txt",header=TRUE,sep='\t')

HVG_allZT_BioVar <- inner_join(AllZT_BioVar,HVG, by=c("Gene"="allZT")) %>%
  select(., Gene, contains("BioVar"))


df <- HVG_allZT_BioVar[2:13]

random_df <- apply(df, 2, 
                   function(x){runif(length(x), min(x), (max(x)))})
random_df <- as.data.frame(random_df)
# Standardize the data sets
df <- scale(df)
random_df <- scale(random_df)

fviz_pca_ind(prcomp(df), title = "PCA - Iris data", 
             habillage = "none",  palette = "jco",
             geom = "point", ggtheme = theme_classic(),
             legend = "bottom")

# Plot the random df
fviz_pca_ind(prcomp(random_df), title = "PCA - Random data", 
             geom = "point", ggtheme = theme_classic())

km.res1 <- kmeans(df, 3)
fviz_cluster(list(data = df, cluster = km.res1$cluster),
             ellipse.type = "norm", geom = "point", stand = FALSE,
             palette = "jco", ggtheme = theme_classic())

fviz_dend(hclust(dist(random_df)), k = 3, k_colors = "jco",  
          as.ggplot = TRUE, show_labels = FALSE)

res <- get_clust_tendency(df, n = nrow(df), graph = FALSE)
res$hopkins_stat

info <- hopkins(counts, n = nrow(counts)-1, byrow = F, header = F)
