## Grid based quantification of fluorescent intensities
## Sample script for in vitro mouse ESCs image dataset

## Set your working directory to the unzipped folder
setwd("PATH_TO_DemoScripts")

## Load required packages
library(Rcpp) # R implementation of C++
library(EBImage) # Image I/O and manipulation
library(mclust) # Gaussian finite mixed model
library(beeswarm) # Beeswarm plot
library(lawstat) # Brunner-Munzel Test

## Initial parameter
# Size of grid in number of pixels
g <- 20

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

######################################################################
#   Preparation of dataset                                           #
#   This is a core procedure of GBIQ                                 #
######################################################################

## Load image set to be analysed
# Read 8-bit grey-scale image files
oct4 <- readImage("./esc/oct4_240.tif")
egfp <- readImage("./esc/egfp_240.tif")
dapi <- readImage("./esc/dapi_240.tif")

# Compute median intensity and IQR of the grids
dapi_med <- GridMedian(dapi, g)
egfp_med <- GridMedian(egfp, g)
oct4_med <- GridMedian(oct4, g)
dapi_iqr <- GridIQR(dapi, g)
egfp_iqr <- GridIQR(egfp, g)
oct4_iqr <- GridIQR(oct4, g)

test <- data.frame(cbind(DAPI_Median=as.vector(dapi_med), DAPI_IQR=as.vector(dapi_iqr), 
                         EGFP_Median=as.vector(egfp_med), EGFP_IQR=as.vector(egfp_iqr), 
                         Oct4_Median=as.vector(oct4_med), Oct4_IQR=as.vector(oct4_iqr)))

## Look at particular grid and plot pixel intensity of the grid to generate Figure 1B
# Extract 10,9-th grid (G10,9)
dapi10_9 <- GridCore(dapi, g, 10, 9)
egfp10_9 <- GridCore(egfp, g, 10, 9)
oct410_9 <- GridCore(oct4, g, 10, 9)

display(dapi10_9)
display(egfp10_9)
display(oct410_9)

# Boxplot and beeswarm of the extracted grid for Figure 1B
boxplot(list(as.vector(dapi10_9), as.vector(egfp10_9), as.vector(oct410_9)), names=c("G10,9(H33342)","G10,9(EGFP)","G10,9(anti-Oct4)"), main="Figure 1B", ylim=c(0,1), ylab="fluorescent intensity/pixel")
beeswarm(list(as.vector(dapi10_9), as.vector(egfp10_9), as.vector(oct410_9)), pch=16, cex=0.75, col=rgb(0,0,0,0.25), add=T)

######################################################################
#   Post-process without "Median_IQR filter" procedure               #
#   This procedure was applied to generate Figure 2A-2D              #
######################################################################

## Classification without Median_IQR filter (Figure 2A-D)
# Initial classification of the dataset
mc.test <- Mclust(test[,c(1,3,5)])
summary(mc.test, parameters=T)
class.test <- mc.test$classification

# Scatter plot for Figure 2A
plot(Oct4_Median~EGFP_Median, data=test[class.test==1,], pch=19, xlim=c(0,0.5), ylim=c(0,0.4), col=rgb(0,0,0,0.25), xlab="EGFP intensity/grid", ylab="anti-Oct4 intensity/grid", main="Figure 2A", bty="n")
par(new=T)
plot(Oct4_Median~EGFP_Median, data=test[class.test==3,], pch=19, xlim=c(0,0.5), ylim=c(0,0.4), col=rgb(0.5,1,0.25,0.5), ann=F, axes=F)
par(new=T)
plot(Oct4_Median~EGFP_Median, data=test[class.test==2,], pch=19, xlim=c(0,0.5), ylim=c(0,0.4), col=rgb(1,0.5,0.5,0.5), ann=F, axes=F)
par(new=T)
plot(Oct4_Median~EGFP_Median, data=test[class.test==4,], pch=19, xlim=c(0,0.5), ylim=c(0,0.4), col=rgb(0,0.5,0.75,0.5), ann=F, axes=F)
legend("bottomright", legend=c("class 1","class 2","class 3","class 4"), pch=19, col=c(rgb(0,0,0,0.5),rgb(1,0.5,0.5,0.5),rgb(0.5,1,0.25,0.5),rgb(0,0.5,0.75,0.5)), bty="n")

# Scatter plot for Figure 2B
plot(Oct4_Median~DAPI_Median, data=test[class.test==1,], pch=19, xlim=c(0,0.7), ylim=c(0,0.4), col=rgb(0,0,0,0.25), xlab="H33342 intensity/grid", ylab="anti-Oct4 intensity/grid", main="Figure 2B", bty="n")
par(new=T)
plot(Oct4_Median~DAPI_Median, data=test[class.test==3,], pch=19, xlim=c(0,0.7), ylim=c(0,0.4), col=rgb(0.5,1,0.25,0.5), ann=F, axes=F)
par(new=T)
plot(Oct4_Median~DAPI_Median, data=test[class.test==2,], pch=19, xlim=c(0,0.7), ylim=c(0,0.4), col=rgb(1,0.5,0.5,0.5), ann=F, axes=F)
par(new=T)
plot(Oct4_Median~DAPI_Median, data=test[class.test==4,], pch=19, xlim=c(0,0.7), ylim=c(0,0.4), col=rgb(0,0.5,0.75,0.5), ann=F, axes=F)
legend("topleft", legend=c("class 1","class 2","class 3","class 4"), pch=19, col=c(rgb(0,0,0,0.5),rgb(1,0.5,0.5,0.5),rgb(0.5,1,0.25,0.5),rgb(0,0.5,0.75,0.5)), bty="n")

# Scatter plot for EGFP~H33342
plot(EGFP_Median~DAPI_Median, data=test[class.test==1,], pch=19, xlim=c(0,0.7), ylim=c(0,0.5), col=rgb(0,0,0,0.25), xlab="H33342 intensity/grid", ylab="EGFP intensity/grid", main="EGFP~H33342", bty="n")
par(new=T)
plot(EGFP_Median~DAPI_Median, data=test[class.test==3,], pch=19, xlim=c(0,0.7), ylim=c(0,0.5), col=rgb(0.5,1,0.25,0.5), ann=F, axes=F)
par(new=T)
plot(EGFP_Median~DAPI_Median, data=test[class.test==2,], pch=19, xlim=c(0,0.7), ylim=c(0,0.5), col=rgb(1,0.5,0.5,0.5), ann=F, axes=F)
par(new=T)
plot(EGFP_Median~DAPI_Median, data=test[class.test==4,], pch=19, xlim=c(0,0.7), ylim=c(0,0.4), col=rgb(0,0.5,0.75,0.5), ann=F, axes=F)
legend("topleft", legend=c("class 1","class 2","class 3","class 4"), pch=19, col=c(rgb(0,0,0,0.5),rgb(1,0.5,0.5,0.5),rgb(0.5,1,0.25,0.5),rgb(0,0.5,0.75,0.5)), bty="n")

# Density plot for Figure 2C
plot(density(test[class.test==1,]$EGFP_Median,bw=0.05),xlim=c(-0.1,0.6),ylim=c(0,8),col=rgb(0.5,0.5,0.5),xlab="EGFP intensity/grid",ylab="kernel density",main="Figure 2C",bty="n")
par(new=T)
plot(density(test[class.test==2,]$EGFP_Median,bw=0.05),xlim=c(-0.1,0.6),ylim=c(0,8),col=rgb(1,0.5,0.5),ann=F,axes=F)
par(new=T)
plot(density(test[class.test==3,]$EGFP_Median,bw=0.05),xlim=c(-0.1,0.6),ylim=c(0,8),col=rgb(0.5,1,0.25),ann=F,axes=F)
par(new=T)
plot(density(test[class.test==4,]$EGFP_Median,bw=0.05),xlim=c(-0.1,0.6),ylim=c(0,8),col=rgb(0,0.5,0.75),ann=F,axes=F)
rug(test[class.test==1,]$EGFP_Median,col=rgb(0.5,0.5,0.5,0.5))
rug(test[class.test==2,]$EGFP_Median,col=rgb(1,0.5,0.5,0.5))
rug(test[class.test==3,]$EGFP_Median,col=rgb(0.5,1,0.25,0.5))
rug(test[class.test==4,]$EGFP_Median,col=rgb(0,0.5,0.75,0.5))
legend("topright", legend=c("class 1","class 2","class 3","class 4"), lty="solid", col=c(rgb(0.5,0.5,0.5),rgb(1,0.5,0.5),rgb(0.5,1,0.25),rgb(0,0.5,0.75)), bty="n")

# Density plot for Figure 2D
plot(density(test[class.test==1,]$Oct4_Median,bw=0.05),xlim=c(-0.1,0.5),ylim=c(0,8),col=rgb(0.5,0.5,0.5),xlab="anti-Oct4 intensity/grid",ylab="kernel density",main="Figure 2D",bty="n")
par(new=T)
plot(density(test[class.test==2,]$Oct4_Median,bw=0.05),xlim=c(-0.1,0.5),ylim=c(0,8),col=rgb(1,0.5,0.5),ann=F,axes=F)
par(new=T)
plot(density(test[class.test==3,]$Oct4_Median,bw=0.05),xlim=c(-0.1,0.5),ylim=c(0,8),col=rgb(0.5,1,0.25),ann=F,axes=F)
par(new=T)
plot(density(test[class.test==4,]$Oct4_Median,bw=0.05),xlim=c(-0.1,0.5),ylim=c(0,8),col=rgb(0,0.5,0.75),ann=F,axes=F)
rug(test[class.test==1,]$Oct4_Median,col=rgb(0.5,0.5,0.5,0.5))
rug(test[class.test==2,]$Oct4_Median,col=rgb(1,0.5,0.5,0.5))
rug(test[class.test==3,]$Oct4_Median,col=rgb(0.5,1,0.25,0.5))
rug(test[class.test==4,]$Oct4_Median,col=rgb(0,0.5,0.75,0.5))
legend("topright", legend=c("class 1","class 2","class 3","class 4"), lty="solid", col=c(rgb(0.5,0.5,0.5),rgb(1,0.5,0.5),rgb(0.5,1,0.25),rgb(0,0.5,0.75)), bty="n")

# Density plot for intensity of H33342
plot(density(test[class.test==1,]$DAPI_Median,bw=0.05),xlim=c(-0.1,0.9),ylim=c(0,8),col=rgb(0.5,0.5,0.5),xlab="H33342 intensity/grid",ylab="kernel density",main="H33342 intensity",bty="n")
par(new=T)
plot(density(test[class.test==2,]$DAPI_Median,bw=0.05),xlim=c(-0.1,0.9),ylim=c(0,8),col=rgb(1,0.5,0.5),ann=F,axes=F)
par(new=T)
plot(density(test[class.test==3,]$DAPI_Median,bw=0.05),xlim=c(-0.1,0.9),ylim=c(0,8),col=rgb(0.5,1,0.25),ann=F,axes=F)
par(new=T)
plot(density(test[class.test==4,]$DAPI_Median,bw=0.05),xlim=c(-0.1,0.9),ylim=c(0,8),col=rgb(0,0.5,0.75),ann=F,axes=F)
rug(test[class.test==1,]$DAPI_Median,col=rgb(0.5,0.5,0.5,0.5))
rug(test[class.test==2,]$DAPI_Median,col=rgb(1,0.5,0.5,0.5))
rug(test[class.test==3,]$DAPI_Median,col=rgb(0.5,1,0.25,0.5))
rug(test[class.test==4,]$DAPI_Median,col=rgb(0,0.5,0.75,0.5))
legend("topright", legend=c("class 1","class 2","class 3","class 4"), lty="solid", col=c(rgb(0.5,0.5,0.5),rgb(1,0.5,0.5),rgb(0.5,1,0.25),rgb(0,0.5,0.75)), bty="n")

######################################################################
#   Introducing "Median_IQR filter"                                  #
#   This procedure was applied to generate Figure 2E & 2F            #
######################################################################

## Classification with Median_IQR filter (Figure 2E & 2F)
# Classify the grids according to H33342 median intensity and its IQR (Median_IQR filter)
mc.med_iqr <- Mclust(test[,1:2])
mean(test$DAPI_Median);mean(test$DAPI_IQR) # Overall mean values of median intensity and IQR
summary(mc.med_iqr, parameters=T)
test <- cbind(test, Median_IQR=mc.med_iqr$classification)
test.class <- which((mc.med_iqr$parameters$mean[1,]>mean(test$DAPI_Median))*(mc.med_iqr$parameters$mean[2,]<mean(test$DAPI_IQR))==1, arr.ind=TRUE)
test.class # Class label for Median_IQR filtered sub-class

# Display the classified grids
mask1 <- as.integer(mc.med_iqr$classification==1)
mask2 <- as.integer(mc.med_iqr$classification==2)
mask3 <- as.integer(mc.med_iqr$classification==3)
mask4 <- as.integer(mc.med_iqr$classification==4)

mc.class1_mask <- matrix(mask1, nrow=nrow(dapi_med), ncol=ncol(dapi_med))
mc.class2_mask <- matrix(mask2, nrow=nrow(dapi_med), ncol=ncol(dapi_med))
mc.class3_mask <- matrix(mask3, nrow=nrow(dapi_med), ncol=ncol(dapi_med))
mc.class4_mask <- matrix(mask4, nrow=nrow(dapi_med), ncol=ncol(dapi_med))

display(mc.class1_mask)
display(mc.class2_mask)
display(mc.class3_mask)
display(mc.class4_mask)

# Pick up representative grids of Median_IQR filtered classes for Supplemental Figure S2D
dapi2_10 <- GridCore(dapi, g, 2, 10)
dapi3_3 <- GridCore(dapi, g, 3, 3)
dapi6_7 <- GridCore(dapi, g, 6, 7)
dapi4_12 <- GridCore(dapi, g, 4, 12)

display(dapi2_10)
display(dapi3_3)
display(dapi6_7)
display(dapi4_12)

# Boxplot and beeswarm of the representative grids
boxplot(list(as.vector(dapi2_10), as.vector(dapi3_3), as.vector(dapi6_7), as.vector(dapi4_12)), names=c("G2,10(class 1)","G3,3(class 2)","G6,7(class 3)","G4,12(class 4)"), main="Supplemental Figure S2D", ylim=c(0,1), ylab="fluorescent intensity/pixel")
beeswarm(list(as.vector(dapi2_10), as.vector(dapi3_3), as.vector(dapi6_7), as.vector(dapi4_12)), pch=16, cex=0.75, col=rgb(0,0,0,0.25), add=T)

## Classiy the filtered grids utilizing H33342, EGFP and anti-Oct4 fluorescent intensities
mc.test <- Mclust(test[test$Median_IQR==test.class,c(1,3,5)])
summary(mc.test, parameters=T)

# Display the classified grids
key <- as.integer(names(mc.test$classification))
mask <- rep(0, nrow(test))
for (i in 1:length(key)) {
  mask[key[i]] <- mc.test$classification[i]
}
test <- cbind(test, Class=mask)

# Scatter plot for Figure 2E
plot(Oct4_Median~EGFP_Median, data=test[test$Class==1,], pch=19, xlim=c(0,0.5), ylim=c(0,0.4), col=rgb(1,0.5,0.5,0.5), xlab="EGFP intensity/grid", ylab="anti-Oct4 intensity/grid", main="Figure 2E", bty="n")
par(new=T)
plot(Oct4_Median~EGFP_Median, data=test[test$Class==2,], pch=19, xlim=c(0,0.5), ylim=c(0,0.4), col=rgb(0.5,1,0.25,0.5), ann=F, axes=F)
legend("bottomright", legend=c("class 3 EGFP(+): 51 grids","class 3 EGFP(-): 16 grids"), pch=19, col=c(rgb(1,0.5,0.5,0.5),rgb(0.5,1,0.25,0.5)), bty="n")

# Density plot for Figure 2F
plot(density(test[test$Class==1,5],bw=0.05),xlim=c(-0.1,0.5),ylim=c(0,8),col=rgb(1,0.5,0.5),xlab="anti-Oct4 intensity/grid",ylab="kernel density",main="Figure 2F",bty="n")
par(new=T)
plot(density(test[test$Class==2,5],bw=0.05),xlim=c(-0.1,0.5),ylim=c(0,8),col=rgb(0.5,1,0.25),ann=F,axes=F)
rug(test[test$Class==1,5],col=rgb(1,0.5,0.5,0.5))
rug(test[test$Class==2,5],col=rgb(0.5,1,0.25,0.5))
legend("topleft", legend=c("class3 EGFP(+)","class3 EGFP(-)"), lty="solid", col=c(rgb(1,0.5,0.5),rgb(0.5,1,0.25)), bty="n")

# Brunner-Munzel tests between EGFP(+) and EGFP(-) sub-classes
brunner.munzel.test(test[test$Class==1,1],test[test$Class==2,1]) # H33342 intensity
brunner.munzel.test(test[test$Class==1,3],test[test$Class==2,3]) # EGFP intensity
brunner.munzel.test(test[test$Class==1,5],test[test$Class==2,5]) # anti-Oct4 intensity

# Kernel density estimation of anti-Oct4 intensities between EGFP(+) and EGFP(-) sub-classes
egfp.pos <- density(test[test$Class==1,5],bw=0.05) # EGFP(+)
egfp.neg <- density(test[test$Class==2,5],bw=0.05) # EGFP(-)
egfp.pos$x[egfp.pos$y==max(egfp.pos$y)] # Estimate anti-Oct4 peak of EGFP(+) sub-class
egfp.neg$x[egfp.neg$y==max(egfp.neg$y)] # Estimate anti-Oct4 peak of EGFP(-) sub-class
egfp.pos$x[egfp.pos$y==max(egfp.pos$y)]/egfp.neg$x[egfp.neg$y==max(egfp.neg$y)] # Difference in anti-Oct4 intensity between EGFP(+) and EGFP(-) sub-classes

## END OF SCRIPT
