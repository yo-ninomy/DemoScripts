## Grid based image quantification (GBIQ) of fluorescent intensities
## Sample script for mouse PSM tissue section image dataset

## Set your working directory to the unzipped folder
setwd("PATH_TO_DemoScripts")

## Load required packages
library(Rcpp) # R implementation of C++
library(EBImage) # Image I/O and manipulation
library(kselection) # K-means selection
library(mclust) # Gaussian finite mixed model
library(mgcv) # Generalized additive models (GAM)
library(plot3D) # 3D plot using OpenGL

## Initial parameter
# Size of grid in number of pixels
g <- 16

## Load basic functions for the analysis
# C++ (Rcpp) implementation of core loop
sourceCpp("GridCore.cpp")

# Compute median intensity of g^2 grid
GridMedian <- function(image, g) {
  gn <- dim(image)%/%g
  med <- matrix(0, nrow=gn[1], ncol=gn[2])
  for(i in 1:gn[1]) {
    for(j in 1:gn[2]) {
      med[i,j] = median(GridCore(image=image,g=g,i=i,j=j))
    }
  }
  return(med)
}

# Compute IQR of the grids
GridIQR <- function(image, g) {
  gn <- dim(image)%/%g
  iqr <- matrix(0, nrow=gn[1], ncol=gn[2])
  for(i in 1:gn[1]) {
    for(j in 1:gn[2]) {
      iqr[i,j] = IQR(as.vector(GridCore(image=image,g=g,i=i,j=j)))
    }
  }
  return(iqr)
}

# Estimate number of clusters by k-means
ClusterNumber <- function(data, iter=100, size=100) {
  k <- NULL
  for (i in 1:iter) {
    kn <- num_clusters_all(kselection(data[sample.int(dim(data)[1], size),]))
    k <- c(k, kn)
  }
  return(table(k))
}

######################################################################
#   Preparation of dataset                                           #
#   This is a core procedure of GBIQ                                 #
######################################################################

## Load image set to be analysed
# Read 8-bit grey-scale image files
dapi <- readImage("./psm/dapi_800.tif")
tbx6 <- readImage("./psm/tbx6_800.tif")
mesp <- readImage("./psm/mesp2_800.tif")
rip2 <- readImage("./psm/ripply2_800.tif")

# Compute median intensity and IQR of the grids
dapi_med <- GridMedian(dapi, g)
tbx6_med <- GridMedian(tbx6, g)
mesp_med <- GridMedian(mesp, g)
rip2_med <- GridMedian(rip2, g)
dapi_iqr <- GridIQR(dapi, g)
tbx6_iqr <- GridIQR(tbx6, g)
mesp_iqr <- GridIQR(mesp, g)
rip2_iqr <- GridIQR(rip2, g)

# Look at the grid_median data
display(dapi_med)
display(tbx6_med)
display(mesp_med)
display(rip2_med)

# Prepare a dataset as dataframe format
test <- data.frame(cbind(DAPI_Median=as.vector(dapi_med), DAPI_IQR=as.vector(dapi_iqr), 
                         Tbx6_Median=as.vector(tbx6_med), Tbx6_IQR=as.vector(tbx6_iqr), 
                         Mesp_Median=as.vector(mesp_med), Mesp_IQR=as.vector(mesp_iqr),
                         Rip2_Median=as.vector(rip2_med), Rip2_IQR=as.vector(rip2_iqr)))

######################################################################
#   Single layer of "Median_IQR filter" procedure                    #
#   This procedure was applied to generate Figure 4                  #
######################################################################

## Apply "Median_IQR filter" on the dataset
# Mclust by median intensity and IQR of DAPI channel (Median_IQR filter)
mc.dapi.med_iqr <- Mclust(test[,1:2], G=1:15)
mean(test$DAPI_Median);mean(test$DAPI_IQR) # Overall mean values of median intensity and IQR
summary(mc.dapi.med_iqr, parameters=T)
test <- cbind(test, Median_IQR=mc.dapi.med_iqr$classification)
test.class <- which((mc.dapi.med_iqr$parameters$mean[1,]>mean(test$DAPI_Median))*(mc.dapi.med_iqr$parameters$mean[2,]<mean(test$DAPI_IQR))==1, arr.ind=TRUE)
test.class

# Display the Median_IQR filtered grids
obj <- d <- NULL
for (i in 1:length(test[,9])) {
  d <- test[i,9]==4|test[i,9]==5
  obj <- c(obj, d)
}
filter_mask <- matrix(obj, nrow=nrow(dapi_med), ncol=ncol(dapi_med))
display(filter_mask)

## Classiy the filtered grids utilizing Tbx6, Mesp2 and Ripply2 expression levels
km <- ClusterNumber(test[test$Median_IQR==4|test$Median_IQR==5,c(3,5,7)], iter=500, size=length(test[test$Median_IQR==4|test$Median_IQR==5,1])/10)
km
mc.test <- Mclust(test[test$Median_IQR==4|test$Median_IQR==5,c(3,5,7)], G=3)
summary(mc.test, parameters=T) # Class 2 indicates high mean value of both Mesp2 and Ripply2 median intensities

# Display the classified grids
key <- as.integer(names(mc.test$classification))
mask <- rep(0, nrow(test))
for (i in 1:length(key)) {
  mask[key[i]] <- mc.test$classification[i]
}
test <- cbind(test, Class=mask)

mask1 <- as.integer(mask==1)
mask2 <- as.integer(mask==2)
mask3 <- as.integer(mask==3)

mc.class1_mask <- matrix(mask1, nrow=nrow(dapi_med), ncol=ncol(dapi_med))
mc.class2_mask <- matrix(mask2, nrow=nrow(dapi_med), ncol=ncol(dapi_med))
mc.class3_mask <- matrix(mask3, nrow=nrow(dapi_med), ncol=ncol(dapi_med))

display(mc.class1_mask)
display(mc.class2_mask)
display(mc.class3_mask)

## Generalized Additive Model (GAM)
# Susetting the dataset by the class
d1 <- test[test$Class==1,c(3,5,7)]
d2 <- test[test$Class==2,c(3,5,7)]
d3 <- test[test$Class==3,c(3,5,7)]
colnames(d1) <- colnames(d2) <- colnames(d3) <- c("Tbx6", "Mesp2", "Ripply2")

# 3D plot of these 3 classes
scatter3D(d2$Mesp2, d2$Ripply2, d2$Tbx6, xlim=c(0,1), ylim=c(0,1), zlim=c(0,1), main="mouse PSM, 800 X 800 pix, Single layer of Median_IQR filter", xlab="Mesp2", ylab="Ripply2", zlab="Tbx6", pch=16, col=rgb(1,0.5,0.5))
scatter3D(d1$Mesp2, d1$Ripply2, d1$Tbx6, xlim=c(0,1), ylim=c(0,1), zlim=c(0,1), pch=16, col=rgb(0,0.5,0.75), add=T, bty="n")
scatter3D(d3$Mesp2, d3$Ripply2, d3$Tbx6, xlim=c(0,1), ylim=c(0,1), zlim=c(0,1), pch=16, col=rgb(0.5,1,0.25), add=T, bty="n")

# GAM analysis of the class 2
d2.gam_tbx6 <- gam(Tbx6~s(Mesp2)+s(Ripply2)+s(Mesp2, Ripply2), data=d2)
vis.gam(d2.gam_tbx6, color="cm", view=c("Mesp2", "Ripply2"), main="mouse PSM, 800 X 800 pix, class 2", sub="Tbx6~s(Mesp2)+s(Ripply2)+s(Mesp2, Ripply2)", theta=-55, phi=45)
summary(d2.gam_tbx6)

######################################################################
#   Two layers of "Median_IQR filter" procedure                      #
#   This procedure was applied to generate Supplementary Figure S6   #
######################################################################

# Prepare a dataset as dataframe format
test2 <- data.frame(cbind(DAPI_Median=as.vector(dapi_med), DAPI_IQR=as.vector(dapi_iqr), 
                          Tbx6_Median=as.vector(tbx6_med), Tbx6_IQR=as.vector(tbx6_iqr), 
                          Mesp_Median=as.vector(mesp_med), Mesp_IQR=as.vector(mesp_iqr),
                          Rip2_Median=as.vector(rip2_med), Rip2_IQR=as.vector(rip2_iqr)))

## Apply "Median_IQR filter" on the dataset
# Apply 1st layer of Median_IQR filter on the dataset
mc.dapi.med_iqr.l1 <- Mclust(test2[,1:2], G=2)
summary(mc.dapi.med_iqr.l1, parameters=T)
# Set class label 0 for void and 1 for cell nucleus
# test2 <- cbind(test2, Median_IQR=mc.dapi.med_iqr.l1$classification-1)
test2 <- cbind(test2, Median_IQR=(mc.dapi.med_iqr.l1$classification-2)*-1)

# Apply 2nd layer of Median_IQR filter on the filtered dataset
mc.dapi.med_iqr.l2 <- Mclust(test2[test2$Median_IQR==1,1:2])
mean(test2[test2$Median_IQR==1,]$DAPI_Median);mean(test2[test2$Median_IQR==1,]$DAPI_IQR)  # Overall mean values of median intensity and IQR
summary(mc.dapi.med_iqr.l2, parameters=T)
test.class <- which((mc.dapi.med_iqr.l2$parameters$mean[1,]>mean(test2[test2$Median_IQR==1,]$DAPI_Median))*(mc.dapi.med_iqr.l2$parameters$mean[2,]<mean(test2[test2$Median_IQR==1,]$DAPI_IQR))==1, arr.ind=TRUE)
test.class # Class label(s) for Median_IQR filtered sub-class(es)

# Add class label to the dataset
key <- as.integer(names(mc.dapi.med_iqr.l2$classification))
mask <- rep(0, nrow(test2))
for (i in 1:length(key)) {
  mask[key[i]] <- mc.dapi.med_iqr.l2$classification[i]
}
test2 <- cbind(test2[,-9], Median_IQR=mask)

# Display the Median_IQR filtered grids
obj <- d <- NULL
for (i in 1:nrow(test2)) {
  d <- test2[i,9]==test.class
  obj <- c(obj, d)
}
filter_mask <- matrix(obj, nrow=nrow(dapi_med), ncol=ncol(dapi_med))
display(filter_mask)

# Classiy the filtered grids utilizing Tbx6, Mesp2 and Ripply2 expression levels
km <- ClusterNumber(test2[test2$Median_IQR==test.class,c(3,5,7)], iter=500, size=length(test2[test2$Median_IQR==test.class,1])/10)
km
mc.test <- Mclust(test2[test2$Median_IQR==test.class,c(3,5,7)], G=3)
summary(mc.test, parameters=T) # Class 2 indicates high mean value of both Mesp2 and Ripply2 median intensities

# Display the classified grids
key <- as.integer(names(mc.test$classification))
mask <- rep(0, length(test2[,9]))
for (i in 1:length(key)) {
  mask[key[i]] <- mc.test$classification[i]
}
test2 <- cbind(test2, Class=mask)

mask1 <- as.integer(mask==1)
mask2 <- as.integer(mask==2)
mask3 <- as.integer(mask==3)

mc.class1_mask <- matrix(mask1, nrow=nrow(dapi_med), ncol=ncol(dapi_med))
mc.class2_mask <- matrix(mask2, nrow=nrow(dapi_med), ncol=ncol(dapi_med))
mc.class3_mask <- matrix(mask3, nrow=nrow(dapi_med), ncol=ncol(dapi_med))

display(mc.class1_mask)
display(mc.class2_mask)
display(mc.class3_mask)

## Generalized Additive Model (GAM)
# Susetting the dataset by the class
d1 <- test2[test2$Class==1,c(3,5,7)]
d2 <- test2[test2$Class==2,c(3,5,7)]
d3 <- test2[test2$Class==3,c(3,5,7)]
colnames(d1) <- colnames(d2) <- colnames(d3) <- c("Tbx6", "Mesp2", "Ripply2")

# 3D plot of these 3 classes
scatter3D(d2$Mesp2, d2$Ripply2, d2$Tbx6, xlim=c(0,1), ylim=c(0,1), zlim=c(0,1), main="mouse PSM, 800 X 800 pix, Dual layer2 of Median_IQR filter", xlab="Mesp2", ylab="Ripply2", zlab="Tbx6", pch=16, col=rgb(1,0.5,0.5))
scatter3D(d1$Mesp2, d1$Ripply2, d1$Tbx6, xlim=c(0,1), ylim=c(0,1), zlim=c(0,1), pch=16, col=rgb(0,0.5,0.75), add=T, bty="n")
scatter3D(d3$Mesp2, d3$Ripply2, d3$Tbx6, xlim=c(0,1), ylim=c(0,1), zlim=c(0,1), pch=16, col=rgb(0.5,1,0.25), add=T, bty="n")

# GAM analysis of class 2
d2.gam_tbx6 <- gam(Tbx6~s(Mesp2)+s(Ripply2)+s(Mesp2, Ripply2), data=d2)
vis.gam(d2.gam_tbx6, color="cm", view=c("Mesp2", "Ripply2"), main="mouse PSM, 800 X 800 pix, class 2", sub="Tbx6~s(Mesp2)+s(Ripply2)+s(Mesp2, Ripply2)", theta=-55, phi=45)
summary(d2.gam_tbx6)

## END OF SCRIPT
