library(tidyverse)
library(clustertend)

AllZT <- read.table("lab/Heatmaps_CV2_expression/BrenneckeMod_allZT_RUV2.txt",
                    header=TRUE, sep='\t')

AllZT_BioVar <- select(AllZT, Gene, contains("BioVar"))

AllZT_BioVar[is.na(AllZT_BioVar)]<-0

HVG <- read.table("lab/Heatmaps_CV2_expression/HVG_allZT_nbZT.txt",header=TRUE,sep='\t')

HVG_allZT_BioVar <- inner_join(AllZT_BioVar,HVG, by=c("Gene"="allZT")) %>%
  select(., Gene, contains("BioVar"))

#----------------------Hopkins------------------------------

x <- HVG_allZT_BioVar[2:13]
AllZT_BioVar <- AllZT_BioVar[2:13]

hopkins(x, n = 100, byrow = F, header = T)
#$H
#[1] 0.2261627
hopkins(x, n = 250, byrow = F, header = T)
#$H
#[1] 0.2146838
hopkins(x, n = 500, byrow = F, header = T)
#$H
#[1] 0.2181909


#Run these just to see
hopkins(AllZT_BioVar, n = 100, byrow = F, header = T)
#$H
#[1] 0.1306029
hopkins(AllZT_BioVar, n = 250, byrow = F, header = T)
#$H
#[1] 0.1324841
hopkins(AllZT_BioVar, n = 500, byrow = F, header = T)
#$H
#[1] 0.1308883