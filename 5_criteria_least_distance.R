
### Identify fluorescence protein based on root mean distance to reference population


## read files

populations <- c("GFP", "YFP", "mKO2", "mCherry")

for (i in populations) {
  assign(paste0(i), read.table(file = paste0("Test images/Measurements/",i,".txt"))) 
}

for (i in populations) {
  
  file <- get(paste0(i))
  rownames(file) <- 1:nrow(file)
  for (j in 1:ncol(file)) {
    colnames(file)[j] <- 434+10*(j-1)
  }
  colnames(file)[ncol(file)] <- "TD"
  assign(paste0(i), file)
}
rm(file, i, j)


## read data on reference pop (those are generated in 4_spectrum_plot.R)

populations <- c("GFP", "YFP", "mKO2", "mCherry")

for (i in populations) {
  assign(paste0(i, "_ctr_sum"), read.table(file = paste0("Results/Spectra/", i, "_ctr_sum.txt"))) 
}


## calculate root mean distance from mean
## sqrt(mean((test-ref)^2))

GFP_RMD <- matrix(nrow = nrow(GFP), ncol = 5)
GFP_RMD <- as.data.frame(GFP_RMD)
colnames(GFP_RMD) <- c("pop", "GFP", "YFP", "mKO2", "mCherry")
GFP_RMD$pop <- "GFP"

for (i in 1:nrow(GFP)) {
  scaled <- GFP[i,1:32]/ GFP[i, which.max(GFP[i,1:32]) ]
  GFP_RMD[i,2] <- sqrt(  ( sum( (scaled-GFP_ctr_sum[,2]/GFP_ctr_sum[ which.max(GFP_ctr_sum[,2]) ,2])^2 ) ) )
  GFP_RMD[i,3] <- sqrt(  ( sum( (scaled-YFP_ctr_sum[,2]/YFP_ctr_sum[ which.max(YFP_ctr_sum[,2]) ,2])^2 ) ) )
  GFP_RMD[i,4] <- sqrt(  ( sum( (scaled-mKO2_ctr_sum[,2]/mKO2_ctr_sum[ which.max(mKO2_ctr_sum[,2]) ,2])^2 ) ) )
  GFP_RMD[i,5] <- sqrt(  ( sum( (scaled-mCherry_ctr_sum[,2]/mCherry_ctr_sum[ which.max(mCherry_ctr_sum[,2]) ,2])^2 ) ) )
}


YFP_RMD <- matrix(nrow = nrow(YFP), ncol = 5)
YFP_RMD <- as.data.frame(YFP_RMD)
colnames(YFP_RMD) <- c("pop", "GFP", "YFP", "mKO2", "mCherry")
YFP_RMD$pop <- "YFP"

for (i in 1:nrow(YFP)) {
  scaled <- YFP[i,1:32]/ YFP[i, which.max(YFP[i,1:32]) ]
  YFP_RMD[i,2] <- sqrt(  ( sum( (scaled-GFP_ctr_sum[,2]/GFP_ctr_sum[ which.max(GFP_ctr_sum[,2]) ,2])^2 ) ) )
  YFP_RMD[i,3] <- sqrt(  ( sum( (scaled-YFP_ctr_sum[,2]/YFP_ctr_sum[ which.max(YFP_ctr_sum[,2]) ,2])^2 ) ) )
  YFP_RMD[i,4] <- sqrt(  ( sum( (scaled-mKO2_ctr_sum[,2]/mKO2_ctr_sum[ which.max(mKO2_ctr_sum[,2]) ,2])^2 ) ) )
  YFP_RMD[i,5] <- sqrt(  ( sum( (scaled-mCherry_ctr_sum[,2]/mCherry_ctr_sum[ which.max(mCherry_ctr_sum[,2]) ,2])^2 ) ) )
}

mKO2_RMD <- matrix(nrow = nrow(mKO2), ncol = 5)
mKO2_RMD <- as.data.frame(mKO2_RMD)
colnames(mKO2_RMD) <- c("pop", "GFP", "YFP", "mKO2", "mCherry")
mKO2_RMD$pop <- "mKO2"

for (i in 1:nrow(mKO2)) {
  scaled <- mKO2[i,1:32]/ mKO2[i, which.max(mKO2[i,1:32]) ]
  mKO2_RMD[i,2] <- sqrt(  ( sum( (scaled-GFP_ctr_sum[,2]/GFP_ctr_sum[ which.max(GFP_ctr_sum[,2]) ,2])^2 ) ) )
  mKO2_RMD[i,3] <- sqrt(  ( sum( (scaled-YFP_ctr_sum[,2]/YFP_ctr_sum[ which.max(YFP_ctr_sum[,2]) ,2])^2 ) ) )
  mKO2_RMD[i,4] <- sqrt(  ( sum( (scaled-mKO2_ctr_sum[,2]/mKO2_ctr_sum[ which.max(mKO2_ctr_sum[,2]) ,2])^2 ) ) )
  mKO2_RMD[i,5] <- sqrt(  ( sum( (scaled-mCherry_ctr_sum[,2]/mCherry_ctr_sum[ which.max(mCherry_ctr_sum[,2]) ,2])^2 ) ) )
}

mCherry_RMD <- matrix(nrow = nrow(mCherry), ncol = 5)
mCherry_RMD <- as.data.frame(mCherry_RMD)
colnames(mCherry_RMD) <- c("pop", "GFP", "YFP", "mKO2", "mCherry")
mCherry_RMD$pop <- "mCherry"

for (i in 1:nrow(mCherry)) {
  scaled <- mCherry[i,1:32]/ mCherry[i, which.max(mCherry[i,1:32]) ]
  mCherry_RMD[i,2] <- sqrt(  ( sum( (scaled-GFP_ctr_sum[,2]/GFP_ctr_sum[ which.max(GFP_ctr_sum[,2]) ,2])^2 ) ) )
  mCherry_RMD[i,3] <- sqrt(  ( sum( (scaled-YFP_ctr_sum[,2]/YFP_ctr_sum[ which.max(YFP_ctr_sum[,2]) ,2])^2 ) ) )
  mCherry_RMD[i,4] <- sqrt(  ( sum( (scaled-mKO2_ctr_sum[,2]/mKO2_ctr_sum[ which.max(mKO2_ctr_sum[,2]) ,2])^2 ) ) )
  mCherry_RMD[i,5] <- sqrt(  ( sum( (scaled-mCherry_ctr_sum[,2]/mCherry_ctr_sum[ which.max(mCherry_ctr_sum[,2]) ,2])^2 ) ) )
}

RMD_mean <- rbind(GFP_RMD, YFP_RMD, mKO2_RMD, mCherry_RMD)





## calculate root mean distance from median

GFP_RMD <- matrix(nrow = nrow(GFP), ncol = 5)
GFP_RMD <- as.data.frame(GFP_RMD)
colnames(GFP_RMD) <- c("pop", "GFP", "YFP", "mKO2", "mCherry")
GFP_RMD$pop <- "GFP"

for (i in 1:nrow(GFP)) {
  scaled <- GFP[i,1:32]/ GFP[i, which.max(GFP[i,1:32]) ]
  GFP_RMD[i,2] <- sqrt(  ( sum( (scaled-GFP_ctr_sum[,5]/GFP_ctr_sum[ which.max(GFP_ctr_sum[,5]) ,5])^2 ) ) )
  GFP_RMD[i,3] <- sqrt(  ( sum( (scaled-YFP_ctr_sum[,5]/YFP_ctr_sum[ which.max(YFP_ctr_sum[,5]) ,5])^2 ) ) )
  GFP_RMD[i,4] <- sqrt(  ( sum( (scaled-mKO2_ctr_sum[,5]/mKO2_ctr_sum[ which.max(mKO2_ctr_sum[,5]) ,5])^2 ) ) )
  GFP_RMD[i,5] <- sqrt(  ( sum( (scaled-mCherry_ctr_sum[,5]/mCherry_ctr_sum[ which.max(mCherry_ctr_sum[,5]) ,5])^2 ) ) )
}


YFP_RMD <- matrix(nrow = nrow(YFP), ncol = 5)
YFP_RMD <- as.data.frame(YFP_RMD)
colnames(YFP_RMD) <- c("pop", "GFP", "YFP", "mKO2", "mCherry")
YFP_RMD$pop <- "YFP"

for (i in 1:nrow(YFP)) {
  scaled <- YFP[i,1:32]/ YFP[i, which.max(YFP[i,1:32]) ]
  YFP_RMD[i,2] <- sqrt(  ( sum( (scaled-GFP_ctr_sum[,5]/GFP_ctr_sum[ which.max(GFP_ctr_sum[,5]) ,5])^2 ) ) )
  YFP_RMD[i,3] <- sqrt(  ( sum( (scaled-YFP_ctr_sum[,5]/YFP_ctr_sum[ which.max(YFP_ctr_sum[,5]) ,5])^2 ) ) )
  YFP_RMD[i,4] <- sqrt(  ( sum( (scaled-mKO2_ctr_sum[,5]/mKO2_ctr_sum[ which.max(mKO2_ctr_sum[,5]) ,5])^2 ) ) )
  YFP_RMD[i,5] <- sqrt(  ( sum( (scaled-mCherry_ctr_sum[,5]/mCherry_ctr_sum[ which.max(mCherry_ctr_sum[,5]) ,5])^2 ) ) )
}

mKO2_RMD <- matrix(nrow = nrow(mKO2), ncol = 5)
mKO2_RMD <- as.data.frame(mKO2_RMD)
colnames(mKO2_RMD) <- c("pop", "GFP", "YFP", "mKO2", "mCherry")
mKO2_RMD$pop <- "mKO2"

for (i in 1:nrow(mKO2)) {
  scaled <- mKO2[i,1:32]/ mKO2[i, which.max(mKO2[i,1:32]) ]
  mKO2_RMD[i,2] <- sqrt(  ( sum( (scaled-GFP_ctr_sum[,5]/GFP_ctr_sum[ which.max(GFP_ctr_sum[,5]) ,5])^2 ) ) )
  mKO2_RMD[i,3] <- sqrt(  ( sum( (scaled-YFP_ctr_sum[,5]/YFP_ctr_sum[ which.max(YFP_ctr_sum[,5]) ,5])^2 ) ) )
  mKO2_RMD[i,4] <- sqrt(  ( sum( (scaled-mKO2_ctr_sum[,5]/mKO2_ctr_sum[ which.max(mKO2_ctr_sum[,5]) ,5])^2 ) ) )
  mKO2_RMD[i,5] <- sqrt(  ( sum( (scaled-mCherry_ctr_sum[,5]/mCherry_ctr_sum[ which.max(mCherry_ctr_sum[,5]) ,5])^2 ) ) )
}

mCherry_RMD <- matrix(nrow = nrow(mCherry), ncol = 5)
mCherry_RMD <- as.data.frame(mCherry_RMD)
colnames(mCherry_RMD) <- c("pop", "GFP", "YFP", "mKO2", "mCherry")
mCherry_RMD$pop <- "mCherry"

for (i in 1:nrow(mCherry)) {
  scaled <- mCherry[i,1:32]/ mCherry[i, which.max(mCherry[i,1:32]) ]
  mCherry_RMD[i,2] <- sqrt(  ( sum( (scaled-GFP_ctr_sum[,5]/GFP_ctr_sum[ which.max(GFP_ctr_sum[,5]) ,5])^2 ) ) )
  mCherry_RMD[i,3] <- sqrt(  ( sum( (scaled-YFP_ctr_sum[,5]/YFP_ctr_sum[ which.max(YFP_ctr_sum[,5]) ,5])^2 ) ) )
  mCherry_RMD[i,4] <- sqrt(  ( sum( (scaled-mKO2_ctr_sum[,5]/mKO2_ctr_sum[ which.max(mKO2_ctr_sum[,5]) ,5])^2 ) ) )
  mCherry_RMD[i,5] <- sqrt(  ( sum( (scaled-mCherry_ctr_sum[,5]/mCherry_ctr_sum[ which.max(mCherry_ctr_sum[,5]) ,5])^2 ) ) )
}

RMD_median <- rbind(GFP_RMD, YFP_RMD, mKO2_RMD, mCherry_RMD)



### count cells

Counts <- data.frame(c(0,0,0,0),
                     c(0,0,0,0),
                     c(0,0,0,0),
                     c(0,0,0,0))
colnames(Counts) <- populations
rownames(Counts) <- populations

for (i in 1:nrow(RMD_mean)) {
  pos <- names(which.min(RMD_mean[i, ] ) )
  Counts[RMD_mean[i,1], pos] <- Counts[RMD_mean[i,1], pos] + 1
}
write.table(Counts, file = "Results/count_RMD_mean.txt", sep = '\t')

Counts <- data.frame(c(0,0,0,0),
                     c(0,0,0,0),
                     c(0,0,0,0),
                     c(0,0,0,0))
colnames(Counts) <- populations
rownames(Counts) <- populations

for (i in 1:nrow(RMD_median)) {
  pos <- names(which.min(RMD_median[i, ] ) )
  Counts[RMD_median[i,1], pos] <- Counts[RMD_median[i,1], pos] + 1
}
write.table(Counts, file = "Results/count_RMD_median.txt", sep = '\t')

































