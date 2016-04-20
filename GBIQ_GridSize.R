## Grid based quantification of fluorescent intensities
## An aid script to get optimal size of grid

## Set your working directory to the unzipped folder
setwd("PATH_TO_DemoScripts")

## Load required packages
library(Rcpp) # R implementation of C++
library(EBImage) # Image I/O and manipulation
library(mclust) # Gaussian finite mixed model

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

## Initial parameters
# Size of grid in number of pix
# Four to six grids covering nucleus would be adequate
# Set candidates of the size as a vector of integers
gs <- c(8,10,12,16,20,24,30,40)

######################################################################
#   Preparation of datasets                                          #
######################################################################

## Load image set to be analysed
# Read 8-bit grey-scale DAPI image files
dapi <- readImage("./esc/dapi_240.tif")

# Following loop calculates GridMedian and GridIQR of each g in the candidate vector 
# and write tables as tab-separated ascii file.
for (i in 1:length(gs)) {
  dapi_med <- GridMedian(dapi, gs[i])
  dapi_iqr <- GridIQR(dapi, gs[i])
  test <- data.frame(cbind(DAPI_Median=as.vector(dapi_med), DAPI_IQR=as.vector(dapi_iqr)))
  mc.test <- Mclust(test)
  test <- cbind(test, Median_IQR=mc.test$classification)
  test.class <- which((mc.test$parameters$mean[1,]>mean(test$DAPI_Median))*(mc.test$parameters$mean[2,]<mean(test$DAPI_IQR))==1, arr.ind=TRUE)
  if (length(test.class)==1) {
    incl <- test$Median_IQR==test.class
  } else if (length(test.class)==2) {
    incl <- test$Median_IQR==test.class[1]|test$Median_IQR==test.class[2]
  } else if (length(test.class)==3) {
    incl <- test$Median_IQR==test.class[1]|test$Median_IQR==test.class[2]|test$Median_IQR==test.class[3]
  } else if (length(test.class)==4) {
    incl <- test$Median_IQR==test.class[1]|test$Median_IQR==test.class[2]|test$Median_IQR==test.class[3]|test$Median_IQR==test.class[4]
  } else if (length(test.class)==5) {
    incl <- test$Median_IQR==test.class[1]|test$Median_IQR==test.class[2]|test$Median_IQR==test.class[3]|test$Median_IQR==test.class[4]|test$Median_IQR==test.class[5]
  } else {
    incl <- rep(FALSE, nrow(test))
  }
  test <- cbind(test, Class=incl)
  write.table(test, paste("test_g", gs[i], ".txt", sep=""), append=F, quote=F, sep="\t", row.names=F)
}

######################################################################
#   For Supplementary Figure S3                                      #
######################################################################

g08 <- read.table("test_g8.txt", header=T, sep="\t", quote="")
g10 <- read.table("test_g10.txt", header=T, sep="\t", quote="")
g12 <- read.table("test_g12.txt", header=T, sep="\t", quote="")
g16 <- read.table("test_g16.txt", header=T, sep="\t", quote="")
g20 <- read.table("test_g20.txt", header=T, sep="\t", quote="")
g24 <- read.table("test_g24.txt", header=T, sep="\t", quote="")
g30 <- read.table("test_g30.txt", header=T, sep="\t", quote="")
g40 <- read.table("test_g40.txt", header=T, sep="\t", quote="")

g08.mc <- Mclust(g08[,1:2], G=max(g08[,3]))
g10.mc <- Mclust(g10[,1:2], G=max(g10[,3]))
g12.mc <- Mclust(g12[,1:2], G=max(g12[,3]))
g16.mc <- Mclust(g16[,1:2], G=max(g16[,3]))
g20.mc <- Mclust(g20[,1:2], G=max(g20[,3]))
g24.mc <- Mclust(g24[,1:2], G=max(g24[,3]))
g30.mc <- Mclust(g30[,1:2], G=max(g30[,3]))
g40.mc <- Mclust(g40[,1:2], G=max(g40[,3]))

table(g08.mc$classification)
g08.mc$parameters$mean
mean(g08$DAPI_Median);mean(g08$DAPI_IQR)
g08.class <- which((g08.mc$parameters$mean[1,]>mean(g08$DAPI_Median))*(g08.mc$parameters$mean[2,]<mean(g08$DAPI_IQR))==1, arr.ind=TRUE)
sum(g08.mc$classification==g08.class)/g08.mc$n

table(g10.mc$classification)
g10.mc$parameters$mean
mean(g10$DAPI_Median);mean(g10$DAPI_IQR)
g10.class <- which((g10.mc$parameters$mean[1,]>mean(g10$DAPI_Median))*(g10.mc$parameters$mean[2,]<mean(g10$DAPI_IQR))==1, arr.ind=TRUE)
sum(g10.mc$classification==g10.class)/g10.mc$n

table(g12.mc$classification)
g12.mc$parameters$mean
mean(g12$DAPI_Median);mean(g12$DAPI_IQR)
g12.class <- which((g12.mc$parameters$mean[1,]>mean(g12$DAPI_Median))*(g12.mc$parameters$mean[2,]<mean(g12$DAPI_IQR))==1, arr.ind=TRUE)
sum(g12.mc$classification==g12.class[1]|g12.mc$classification==g12.class[2])/g12.mc$n
# Two out of 6 classes fulfill the criteria. Make a subset of these classes together.
g12.class4_5 <- g12[g12.mc$classification==g12.class[1]|g12.mc$classification==g12.class[2],]

table(g16.mc$classification)
g16.mc$parameters$mean
mean(g16$DAPI_Median);mean(g16$DAPI_IQR)
g16.class <- which((g16.mc$parameters$mean[1,]>mean(g16$DAPI_Median))*(g16.mc$parameters$mean[2,]<mean(g16$DAPI_IQR))==1, arr.ind=TRUE)
sum(g16.mc$classification==g16.class)/g16.mc$n

table(g20.mc$classification)
g20.mc$parameters$mean
mean(g20$DAPI_Median);mean(g20$DAPI_IQR)
g20.class <- which((g20.mc$parameters$mean[1,]>mean(g20$DAPI_Median))*(g20.mc$parameters$mean[2,]<mean(g20$DAPI_IQR))==1, arr.ind=TRUE)
sum(g20.mc$classification==g20.class)/g20.mc$n

table(g24.mc$classification)
g24.mc$parameters$mean
mean(g24$DAPI_Median);mean(g24$DAPI_IQR)
g24.class <- which((g24.mc$parameters$mean[1,]>mean(g24$DAPI_Median))*(g24.mc$parameters$mean[2,]<mean(g24$DAPI_IQR))==1, arr.ind=TRUE)
sum(g24.mc$classification==g24.class)/g24.mc$n

table(g30.mc$classification)
g30.mc$parameters$mean
mean(g30$DAPI_Median);mean(g30$DAPI_IQR)
# g30.class <- which((g30.mc$parameters$mean[1,]>mean(g30$DAPI_Median))*(g30.mc$parameters$mean[2,]<mean(g30$DAPI_IQR))==1, arr.ind=TRUE)
# Non of the classes fulfill the criteria. Select a class indicating its mean of median intensities is above the overall mean.
g30.class <- which((g30.mc$parameters$mean[1,]>mean(g30$DAPI_Median))==1, arr.ind=TRUE)
sum(g30.mc$classification==g30.class)/g30.mc$n

table(g40.mc$classification)
g40.mc$parameters$mean
mean(g40$DAPI_Median);mean(g40$DAPI_IQR)
g40.class <- which((g40.mc$parameters$mean[1,]>mean(g40$DAPI_Median))*(g40.mc$parameters$mean[2,]<mean(g40$DAPI_IQR))==1, arr.ind=TRUE)
sum(g40.mc$classification==g40.class)/g40.mc$n

DAPI_Median <- c(mean(g08$DAPI_Median), mean(g10$DAPI_Median), mean(g12$DAPI_Median), mean(g16$DAPI_Median), mean(g20$DAPI_Median), mean(g24$DAPI_Median), mean(g30$DAPI_Median), mean(g40$DAPI_Median))
DAPI_Median_SD <- c(sd(g08$DAPI_Median), sd(g10$DAPI_Median), sd(g12$DAPI_Median), sd(g16$DAPI_Median), sd(g20$DAPI_Median), sd(g24$DAPI_Median), sd(g30$DAPI_Median), sd(g40$DAPI_Median))
DAPI_IQR <- c(mean(g08$DAPI_IQR), mean(g10$DAPI_IQR), mean(g12$DAPI_IQR), mean(g16$DAPI_IQR), mean(g20$DAPI_IQR), mean(g24$DAPI_IQR), mean(g30$DAPI_IQR), mean(g40$DAPI_IQR))
DAPI_IQR_SD <- c(sd(g08$DAPI_IQR), sd(g10$DAPI_IQR), sd(g12$DAPI_IQR), sd(g16$DAPI_IQR), sd(g20$DAPI_IQR), sd(g24$DAPI_IQR), sd(g30$DAPI_IQR), sd(g40$DAPI_IQR))
Class_Number <- c(g08.mc$G, g10.mc$G, g12.mc$G, g16.mc$G, g20.mc$G, g24.mc$G, g30.mc$G, g40.mc$G)
Obs_Number <- c(g08.mc$n, g10.mc$n, g12.mc$n, g16.mc$n, g20.mc$n, g24.mc$n, g30.mc$n, g40.mc$n)
Class_Obs <- c(sum(g08.mc$classification==g08.class), sum(g10.mc$classification==g10.class), sum(g12.mc$classification==g12.class[1]|g12.mc$classification==g12.class[2]), sum(g16.mc$classification==g16.class), sum(g20.mc$classification==g20.class), sum(g24.mc$classification==g24.class), sum(g30.mc$classification==g30.class), sum(g40.mc$classification==g40.class))
Class_Ratio <- c(sum(g08.mc$classification==g08.class)/g08.mc$n, sum(g10.mc$classification==g10.class)/g10.mc$n, sum(g12.mc$classification==g12.class[1]|g12.mc$classification==g12.class[2])/g12.mc$n, sum(g16.mc$classification==g16.class)/g16.mc$n, sum(g20.mc$classification==g20.class)/g20.mc$n, sum(g24.mc$classification==g24.class)/g24.mc$n, sum(g30.mc$classification==g30.class)/g30.mc$n, sum(g40.mc$classification==g40.class)/g40.mc$n)
Class_Median <- c(g08.mc$parameters$mean[1,g08.class], g10.mc$parameters$mean[1,g10.class], mean(g12.class4_5[,1]), g16.mc$parameters$mean[1,g16.class], g20.mc$parameters$mean[1,g20.class], g24.mc$parameters$mean[1,g24.class], g30.mc$parameters$mean[1,g30.class], g40.mc$parameters$mean[1,g40.class])
Class_Median_SD <- c(sd(g08[g08.mc$classification==g08.class,1]), sd(g10[g10.mc$classification==g10.class,1]), sd(g12.class4_5[,1]), sd(g16[g16.mc$classification==g16.class,1]), sd(g20[g20.mc$classification==g20.class,1]), sd(g24[g24.mc$classification==g24.class,1]), sd(g30[g30.mc$classification==g30.class,1]), sd(g40[g40.mc$classification==g40.class,1]))
Class_IQR <- c(g08.mc$parameters$mean[2,g08.class], g10.mc$parameters$mean[2,g10.class], mean(g12.class4_5[,2]), g16.mc$parameters$mean[2,g16.class], g20.mc$parameters$mean[2,g20.class], g24.mc$parameters$mean[2,g24.class], g30.mc$parameters$mean[2,g30.class], g40.mc$parameters$mean[2,g40.class])
Class_IQR_SD <- c(sd(g08[g08.mc$classification==g08.class,2]), sd(g10[g10.mc$classification==g10.class,2]), sd(g12.class4_5[,2]), sd(g16[g16.mc$classification==g16.class,2]), sd(g20[g20.mc$classification==g20.class,2]), sd(g24[g24.mc$classification==g24.class,2]), sd(g30[g30.mc$classification==g30.class,2]), sd(g40[g40.mc$classification==g40.class,2]))
Eval_Matrix <- as.data.frame(cbind(DAPI_Median, DAPI_Median_SD, DAPI_IQR, DAPI_IQR_SD, Obs_Number, Class_Number, Class_Obs, Class_Ratio, Class_Median, Class_Median_SD, Class_IQR, Class_IQR_SD))
rownames(Eval_Matrix) <- c("g08", "g10", "g12", "g16", "g20", "g24", "g30", "g40")

######################################################################
#   Draw Supplementary Figure S3H                                    #
######################################################################

xaxis <- 1:nrow(Eval_Matrix)
plot(Eval_Matrix$DAPI_Median, ylim=c(0,0.6), type="b", pch=16, col=rgb(0,0.5,0.75), xlab="grid size", ylab="H33342 intensity", main="Effects of grid size", bty="n")
arrows(xaxis, Eval_Matrix$DAPI_Median, xaxis, Eval_Matrix$DAPI_Median+Eval_Matrix$DAPI_Median_SD, angle=90, length=0.1, col=rgb(0,0.5,0.75,0.5))
arrows(xaxis, Eval_Matrix$DAPI_Median, xaxis, Eval_Matrix$DAPI_Median-Eval_Matrix$DAPI_Median_SD, angle=90, length=0.1, col=rgb(0,0.5,0.75,0.5))
par(new=T)
plot(Eval_Matrix$DAPI_IQR, ylim=c(0,0.6), type="b", pch=15, col=rgb(0,0.5,0.75), ann=F, axes=F)
arrows(xaxis, Eval_Matrix$DAPI_IQR, xaxis, Eval_Matrix$DAPI_IQR+Eval_Matrix$DAPI_IQR_SD, angle=90, length=0.1, col=rgb(0,0.5,0.75,0.5))
arrows(xaxis, Eval_Matrix$DAPI_IQR, xaxis, Eval_Matrix$DAPI_IQR-Eval_Matrix$DAPI_IQR_SD, angle=90, length=0.1, col=rgb(0,0.5,0.75,0.5))
par(new=T)
plot(Eval_Matrix$Class_Median, ylim=c(0,0.6), type="b", pch=16, col=rgb(1,0.5,0.5), ann=F, axes=F)
arrows(xaxis, Eval_Matrix$Class_Median, xaxis, Eval_Matrix$Class_Median+Eval_Matrix$Class_Median_SD, angle=90, length=0.1, col=rgb(1,0.5,0.5,0.5))
arrows(xaxis, Eval_Matrix$Class_Median, xaxis, Eval_Matrix$Class_Median-Eval_Matrix$Class_Median_SD, angle=90, length=0.1, col=rgb(1,0.5,0.5,0.5))
par(new=T)
plot(Eval_Matrix$Class_IQR, ylim=c(0,0.6), type="b", pch=15, col=rgb(1,0.5,0.5), ann=F, axes=F)
arrows(xaxis, Eval_Matrix$Class_IQR, xaxis, Eval_Matrix$Class_IQR+Eval_Matrix$Class_IQR_SD, angle=90, length=0.1, col=rgb(1,0.5,0.5,0.5))
arrows(xaxis, Eval_Matrix$Class_IQR, xaxis, Eval_Matrix$Class_IQR-Eval_Matrix$Class_IQR_SD, angle=90, length=0.1, col=rgb(1,0.5,0.5,0.5))
legend("bottomright", legend=c("H33342 median intensity - all grids","IQR of H33342 intensity - all grids","H33342 median intensity - Selected class","IQR of H33342 intensity - Selected class"), pch=c(16,15,16,15), col=c(rgb(0,0.5,0.75),rgb(0,0.5,0.75),rgb(1,0.5,0.5),rgb(1,0.5,0.5)), bty="n")

######################################################################
#   Draw Supplementary Figure S3I                                    #
######################################################################

plot(log2(Eval_Matrix$Obs_Number), ylim=c(4,10), type="b", pch=16, col=rgb(0,0.5,0.75), xlab="grid size", ylab="number of grids (log2)", main="Effects of grid size", bty="n")
par(new=T)
plot(log2(Eval_Matrix$Class_Obs), ylim=c(4,10), type="b", pch=16, col=rgb(1,0.5,0.5), ann=F, axes=F)
abline(h=log2(41), lty="solid")
abline(h=log2(82), lty="dotted")
legend("topright", legend=c("Number of all grids","Number of grids in selected class"), pch=16, col=c(rgb(0,0.5,0.75),rgb(1,0.5,0.5)), bty="n")

######################################################################
#   A matrix of descriptive statistics for each selected class       #
######################################################################

Result <- Eval_Matrix[,c(1,3,5,7,9,11)]
Result
#     DAPI_Median   DAPI_IQR Obs_Number Class_Obs Class_Median  Class_IQR
# g08   0.3301111 0.09922549        900       230    0.4345268 0.05641323
# g10   0.3292858 0.11452546        576       176    0.4377934 0.07482120
# g12   0.3289216 0.12341176        400       144    0.4322032 0.09832516
# g16   0.3295599 0.13777342        225       110    0.4328190 0.11139259
# g20   0.3305964 0.15990605        144        67    0.4547339 0.13321143
# g24   0.3292549 0.17129412        100        55    0.4502235 0.16122725
# g30   0.3334252 0.18455882         64        43    0.4334760 0.19883792
# g40   0.3382353 0.19896514         36        21    0.4408494 0.19801131

## END OF SCRIPT
