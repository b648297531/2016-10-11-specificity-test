###########################################################################################

library(ggplot2)
library(reshape2)

populations <- c("CFP", "GFP", "YFP", "mKO2", "mCherry")

## read data
for (i in populations) {
  assign(i, read.table(file = paste0("Test images (individual)/Measurements/",i,".txt"))) ## test images
  assign(paste0(i,"_ctr"), read.table(file = paste0("unmixed/Measurements/",i,".txt"))) ## ctr images
}

## Make compensation table
## 5x5 data.frame 
compensation <- matrix(nrow = 5, ncol = 5)
compensation <- as.data.frame(compensation)
colnames(compensation) <- populations
rownames(compensation) <- populations

## Compute compensation values. 
## for each populations, compensate so that median of the negative channels are ~50
## In a CFP pop, signal in GFP channel is some multiple of CFP signal plus background and so on.
## median(GFP) = k*median(CFP) + background, solve for k

## background varies from channel to channel
background <- c(90, 100, 60, 10, 18)
names(background) <- populations

## read measurements of background signals in each channels
## ** I havent measure those yet!
## ** maybe same formats as other files aka. CFP_ctr_background, GFP_ctr_background, .... 

#for (i in populations) {
#  assign(paste0(i,"_ctr_background"), read.table(file = paste0("unmixed/Measurements/",i,".txt"))) ## ctr images
#}

for (i in populations) {
  ## population to be compensated
  target <- paste0(i, "_ctr")
  pos_median <- median( get(target)$i )
  print(pos_median)
  ## go through each channel and calculate k
  for (j in populations) {
    if (j != i) {
      neg_median <- median( get(target)$j )
      k <- (neg_median-median(get(paste0(i, "_ctr_background"))$j))/pos_median
      compensation[j,i] <- k
    }
    else {
      compensation[j,i] <- 1
    }
  }
}
rm(i,j,k, target,pos_median,neg_median)
write.table(compensation, file = "Results/compensation.txt", sep = '\t')

### Copy data to colour.comp and apply compensations
{
  CFP.comp <- CFP
  for (i in 1:nrow(CFP.comp)) {
    CFP.comp[i,1:5] <- solve(compensation,CFP.comp[i,1:5])
  }
  GFP.comp <- GFP
  for (i in 1:nrow(GFP.comp)) {
    GFP.comp[i,1:5] <- solve(compensation,GFP.comp[i,1:5])
  }
  YFP.comp <- YFP
  for (i in 1:nrow(YFP.comp)) {
    YFP.comp[i,1:5] <- solve(compensation,YFP.comp[i,1:5])
  }
  mKO2.comp <- mKO2
  for (i in 1:nrow(mKO2.comp)) {
    mKO2.comp[i,1:5] <- solve(compensation,mKO2.comp[i,1:5])
  }
  mCherry.comp <- mCherry
  for (i in 1:nrow(mCherry.comp)) {
    mCherry.comp[i,1:5] <- solve(compensation,mCherry.comp[i,1:5])
  }
}

## graph intensity distribution of all channels in all populations
## couldnt figure out how to save them, will save manually
for (i in populations) {
  x <- melt(get(paste0(i,".comp")))
  assign(paste0(i,".comp_boxplot"), ggplot(data = x, aes(x = variable,y = value, fill=variable), legend = FALSE) + geom_boxplot(show.legend = FALSE))
}
rm(x)

## Saving compensated boxplots manually
pdf(file = paste0("Results/", "CFP.comp", "intensity_all_channels_box_median.pdf" ))
CFP.comp_boxplot
dev.off()
pdf(file = paste0("Results/", "GFP.comp", "intensity_all_channels_box_median.pdf" ))
GFP.comp_boxplot
dev.off()
pdf(file = paste0("Results/", "YFP.comp", "intensity_all_channels_box_median.pdf" ))
YFP.comp_boxplot
dev.off()
pdf(file = paste0("Results/", "mKO2.comp", "intensity_all_channels_box_median.pdf" ))
mKO2.comp_boxplot
dev.off()
pdf(file = paste0("Results/", "mCherry.comp", "intensity_all_channels_box_median.pdf" ))
mCherry.comp_boxplot
dev.off()


## Identify each cells as positive or negative in each channels for each populations
## Use lowest 1% of the population as cut off?

### compute cut off point
### 1% false neg rate

Threadsholds <- numeric(length = 5)
names(Threadsholds) <- populations

for (i in 1:5) {
  pop <- paste0(populations[i], ".comp")
  if (round(nrow(get(pop))*0.01, digits = 0) != 0){
    which <- order(get(pop)[,populations[i]], decreasing = FALSE)[round(nrow(get(pop))*0.01, digits = 0)]
    Threadsholds[i] <- get(pop)[which,populations[i]]
  }
  else {
    which <- order(get(pop)[,populations[i]], decreasing = FALSE)[1]
    Threadsholds[i] <- get(pop)[which,populations[i]]
  }
}
rm(pop, i, which)

### go through each events in all popuations and make a call
Counts <- matrix(nrow = 5, ncol = 5)
Counts <- as.data.frame(Counts)
colnames(Counts) <- populations
rownames(Counts) <- populations

for (i in populations) {
  pop <- paste0(i, ".comp")
  number <- c(0,0,0,0,0)
  
  for (j in 1:5) {
    for (k in 1:nrow(get(pop))) {
      if (get(pop)[k,j] >= Threadsholds[i]){
        number[j] <- number[j] +1
      }
    }
  }
  Counts[i,] <- signif(number/nrow(get(pop)), digits = 2)
}
rm(number,i,j,k,pop)

write.table(Counts, file = "Results/count_1%_FN.txt", sep = '\t')

### count again with threadshold of 200
Counts <- matrix(nrow = 5, ncol = 5)
Counts <- as.data.frame(Counts)
colnames(Counts) <- populations
rownames(Counts) <- populations

for (i in populations) {
  pop <- paste0(i, ".comp")
  number <- c(0,0,0,0,0)
  for (j in 1:5) {
    for (k in 1:nrow(get(pop))) {
      if (get(pop)[k,j] >= 500){
        number[j] <- number[j] +1
      }
    }
  }
  Counts[i,] <- round(number/nrow(get(pop)), digits = 2)
}
rm(number,i,j,k,pop)
write.table(Counts, file = "Results/count_500_cutoff.txt", sep = '\t')
