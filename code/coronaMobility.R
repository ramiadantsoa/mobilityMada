#### building the data matrix from the google maps distance ####

# import data
mobG <- read.csv("/Users/tanjona/Box Sync/projects/corona_mada/data/Mobility_GoogleMaps_Matrix_TR.csv", header = TRUE, row.names = 1)

# fill the upper diagonal
for(i in 1:22){
  for(j in (i+1):22){
    mobG[i,j] <- mobG[j,i]
  }
}
mobG <- mobG[, -23]

# convert to seconds
mobG1 <- mobG
for(i in 1:22){
  for(j in 1:22){
    mobG1[i,j]
    if(is.na(mobG1[i,j]) == TRUE | mobG1[i,j] == "0"){
    mobG1[i,j] <- mobG1[j,i]
    }else{
      vec <- as.numeric(strsplit(mobG1[i,j], ":")[[1]])
      newval <- vec[1:2] %*% c(3600, 60)
      mobG1[i,j] <- as.numeric(newval)
    }
  }
}


# what to do with Maintirano
# Tsiroanomandidy to Ambodimanara is 3:35 that is 12900
# And then add two days from Ambodimanara 48*3600 = 172800
# the total is 185700
mobG1$Maintirano <- 185700
mobG1[10,] <- 185700
mobG1[10,10] <- 0

mobG1$Maintirano[1]

# add time to Tsiroanimandidy ID = 6
for(i in 1:22){
  mobG1$Maintirano[i] <- mobG1$Maintirano[i] + as.numeric(mobG1$Tsiroanomandidy[i]) 
}  
mobG1[10,10] <- 0

capital_name <- colnames(mobG1)

reg_name <-  read.csv("/Users/tanjona/Box Sync/projects/corona_mada/data/key_reg_capital.csv", header = FALSE)[,8]

colnames(mobG1) <- reg_name
rownames(mobG1) <- reg_name

#replace the order to match mathematica code
reg_name <- read.csv('/Users/tanjona/Box Sync/projects/corona_mada/data/MDG_5yrs_InternalMigFlows_2010/keys_ID.csv', header = T)

temp <- mobG1

colnames(mobG1) <- reg_name$Mathematica
rownames(mobG1) <- reg_name$Mathematica

# reorder the matrix
for (i in 1:22) {
  for (j in 1:22) {
    xx <- reg_name[i,'Mathematica']
    yy <- reg_name[j,'Mathematica']
    mobG1[i,j] <- temp[xx,yy]
  }
}

mobG2 <- lapply(mobG1, function(x) as.numeric(as.character(x)))


write.csv(mobG2, "/Users/tanjona/Box Sync/projects/corona_mada/data/mobility_google.csv", row.names = TRUE)

#### scaling the mobility matrices ####
mGo <- read.csv("/Users/tanjona/Box Sync/projects/corona_mada/data/mobility_google.csv", row.names = 1)
rownames(mGo) <- reg_name$Mathematica

mGr <- read.csv("/Users/tanjona/Box Sync/projects/corona_mada/data/mobility_gravity.csv", row.names = 1)

mIm <- read.csv("/Users/tanjona/Box Sync/projects/corona_mada/data/mobility_imf.csv", row.names = 1)

mGo <- 1/mGo
diag(mGo) <- 0

mGo <- mGo/(sum(mGo))
mGr <- mGr/(sum(mGr))
mIm <- mIm/(sum(mIm))

mean(as.dist(mGo))
sd(as.dist(mGo))
sd(as.dist(mGo))/mean(as.dist(mGo))

mean(as.dist(mGr))
sd(as.dist(mGr))
sd(as.dist(mGr))/mean(as.dist(mGr))


mean(as.dist(mIm))
sd(as.dist(mIm))
sd(as.dist(mIm))/mean(as.dist(mIm))

library(tidyr)
library(ggplot2)

df <- data.frame(type = rep(c("Google", "Gravity", "InterFlow"), each = 231), dist = c(as.dist(mGo), as.dist(mGr), as.dist(mIm)))

ggplot(df, aes(x = dist))+
  geom_histogram(data=subset(df, type == 'Google'), bins = 20, fill = "red", alpha = 0.5)+
  geom_histogram(data=subset(df, type == 'Gravity'), bins = 20, fill = "black", alpha = 0.5)+
  geom_histogram(data=subset(df, type == 'InterFlow'), bins = 20, fill = "cyan", alpha = 0.5)+
  theme_bw() +
  scale_x_log10()+
  xlab("Distribution of distances")

library(vegan)
library(ggplot2)
require(gclus)

# for google data
clust_mGo <- hclust(as.dist(mGo),method='centroid')
clust_mGo <- reorder.hclust(clust_mGo, mGo)
dmGo <- as.dendrogram(clust_mGo)
heatmap(as.matrix(mGo),Rowv=dmGo,symm=T,col = hcl.colors(20, "Blue-Red 3", rev = TRUE)) 

# for gravity data
clust_mGr <- hclust(as.dist(mGr),method='centroid')
clust_mGr <- reorder.hclust(clust_mGr, mGr)
dmGr <- as.dendrogram(clust_mGr)
heatmap(as.matrix(mGr),Rowv=dmGr,symm=T,col = hcl.colors(20, "Blue-Red 3", rev = TRUE))

# for internal flow data
clust_mIm <- hclust(as.dist(mIm),method='centroid')
clust_mIm <- reorder.hclust(clust_mIm, mIm)
dmIm <- as.dendrogram(clust_mIm)
heatmap(as.matrix(mIm),Rowv=dmIm,symm=T,col = hcl.colors(20, "Blue-Red 3", rev = FALSE))



#### scaling the mobility by column ####
mGo <- read.csv("/Users/tanjona/Box Sync/projects/corona_mada/data/mobility_google.csv", row.names = 1)
rownames(mGo) <- reg_name$Mathematica

mGr <- read.csv("/Users/tanjona/Box Sync/projects/corona_mada/data/mobility_gravity.csv", row.names = 1)

mIm <- read.csv("/Users/tanjona/Box Sync/projects/corona_mada/data/mobility_imf.csv", row.names = 1)

mGo <- 1/mGo
diag(mGo) <- 0
mGo <- mGo/lapply(mGo, sum)
mGr <- mGr/lapply(mGr, sum)
mIm <- mIm/lapply(mIm, sum)


sumGo1 <- round(mean(as.dist(mGo)),3)
sd(as.dist(mGo))
sumGo2 <- round(sd(as.dist(mGo))/mean(as.dist(mGo)),3)

sumGr1 <- round(mean(as.dist(mGr)),3)
sd(as.dist(mGr))
sumGr2 <- round(sd(as.dist(mGr))/mean(as.dist(mGr)),3)


sumIm1 <- round(mean(as.dist(mIm)),3)
sd(as.dist(mIm))
sumIm2 <- round(sd(as.dist(mIm))/mean(as.dist(mIm)),3)

df.sum <- data.frame(mean = c(sumGo1, sumIm1,sumGr1), CV = c(sumGo2, sumIm2, sumGr2))
rownames(df.sum) <- c("Google","InternalMF","Gravity")

library(tidyr)
library(ggplot2)

df <- data.frame(type = rep(c("Google", "Gravity", "InterFlow"), each = 231), dist = c(as.dist(mGo), as.dist(mGr), as.dist(mIm)))

ggplot(df, aes(x = dist))+
  geom_histogram(data=subset(df, type == 'Google'), bins = 15, fill = "red", alpha = 0.8)+
  geom_histogram(data=subset(df, type == 'Gravity'), bins = 15, fill = "black", alpha = 0.4)+
  geom_histogram(data=subset(df, type == 'InterFlow'), bins = 15, fill = "cyan", alpha = 0.5)+
  theme_bw() +
  scale_x_log10()+
  xlab("Distribution of distances")

library(vegan)
library(ggplot2)
require(gclus)


# for google data
clust_mGo <- hclust(as.dist(mGo),method='average')
clust_mGo <- reorder.hclust(clust_mGo, mGo)
dmGo <- as.dendrogram(clust_mGo)
heatmap(as.matrix(mGo),Rowv=dmGo,symm=T,  col = hcl.colors(20, "gray", rev = TRUE), margins = c(10,10)) 

# for gravity data
clust_mGr <- hclust(as.dist(mGr),method='average')
clust_mGr <- reorder.hclust(clust_mGr, mGr)
dmGr <- as.dendrogram(clust_mGr)
heatmap(as.matrix(mGr), Rowv = dmGr,symm = T,col = hcl.colors(20, "gray", rev = TRUE), margins = c(10,10))

# for internal flow data
clust_mIm <- hclust(as.dist(mIm),method='average')
clust_mIm <- reorder.hclust(clust_mIm, mIm)
dmIm <- as.dendrogram(clust_mIm)
heatmap(as.matrix(mIm), Rowv = dmIm,symm = T,col = hcl.colors(20, "gray", rev = TRUE), margins = c(10,10))

