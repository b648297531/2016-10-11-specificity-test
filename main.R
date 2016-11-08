## organize measurement files so that each colours are in one file
## Change arrangement as well
{
populations <- c("CFP", "GFP", "YFP", "mKO2", "mCherry")

### read all the files. 
for (i in populations) {
  for (j in 1:5){
    file_name <- paste0(i, j, ".xls")
    file <- t(read.table(file = paste0("Test images/Measurements/", file_name), 
                       sep = '\t', header = TRUE, row.names = 1))
    rownames(file) <- 1:nrow(file)
    colnames(file) <- c("CFP", "GFP", "YFP", "mKO2", "mCherry", "TD")
    assign(paste0(i, j), value = as.data.frame(file))
  }
}
rm(i,j,file, file_name)

## Combine files
for (i in populations) {
  assign(i, 
         value = rbind(get(paste0(i,1)), get(paste0(i,2)), get(paste0(i,3)), 
                       get(paste0(i,4)), get(paste0(i,5))
                       )
         )
}
## Write out combined files
write.table(CFP, file = "Test images/CFP.txt", sep = '\t')
write.table(GFP, file = "Test images/GFP.txt", sep = '\t')
write.table(YFP, file = "Test images/YFP.txt", sep = '\t')
write.table(mKO2, file = "Test images/mKO2.txt", sep = '\t')
write.table(mCherry, file = "Test images/mCherry.txt", sep = '\t')

## Import measurements on ctrs and transform them to the corrected order.
## save files. 
populations <- c("CFP", "GFP", "YFP", "mKO2", "mCherry")

for (i in populations) {
  file_name <- paste0(i, " ctr.xls")
  file <- t(read.table(file = paste0("unmixed/Measurements/", file_name), 
                       sep = '\t', header = TRUE, row.names = 1))
  rownames(file) <- 1:nrow(file)
  colnames(file) <- c("CFP", "GFP", "YFP", "mKO2", "mCherry", "TD")
  assign(paste0(i, "_ctr"), value = as.data.frame(file))
  
}
rm(i,file_name, file)
## retain only mean measurement (every 3 rows from 1 to end)
GFP_ctr <- GFP_ctr[seq(1,nrow(GFP_ctr),3),]
CFP_ctr <- CFP_ctr[seq(1,nrow(CFP_ctr),3),]
YFP_ctr <- YFP_ctr[seq(1,nrow(YFP_ctr),3),]
mKO2_ctr <- mKO2_ctr[seq(1,nrow(mKO2_ctr),3),]
mCherry_ctr <- mCherry_ctr[seq(1,nrow(mCherry_ctr),3),]

## save files 
write.table(CFP_ctr, file = "unmixed/Measurements/CFP.txt", sep = '\t')
write.table(GFP_ctr, file = "unmixed/Measurements/GFP.txt", sep = '\t')
write.table(YFP_ctr, file = "unmixed/Measurements/YFP.txt", sep = '\t')
write.table(mKO2_ctr, file = "unmixed/Measurements/mKO2.txt", sep = '\t')
write.table(mCherry_ctr, file = "unmixed/Measurements/mCherry.txt", sep = '\t')
}
## Do not need to run those anymore ^

###########################################################################################
###########################################################################################

library(ggplot2)
library(reshape2)

populations <- c("CFP", "GFP", "YFP", "mKO2", "mCherry")

## read data
for (i in populations) {
  assign(i, read.table(file = paste0("Test images (individual)/Measurements/",i,".txt"))) ## test images
  assign(paste0(i,"_ctr"), read.table(file = paste0("unmixed/Measurements/",i,".txt"))) ## ctr images
}

## Loop through and make 2D scatter plots for ctr
{
for (i in populations) {
  pdf(file = paste0("Results/ctr 2D scatter/", i ,"_positive_", "_ctr.pdf"))
  par(mfrow = c(2,2))
  pop <- populations[populations != i]
  plot(get(paste0(i, "_ctr"))[,c(i, pop[1])])
  plot(get(paste0(i, "_ctr"))[,c(i, pop[2])])
  plot(get(paste0(i, "_ctr"))[,c(i, pop[3])])
  plot(get(paste0(i, "_ctr"))[,c(i, pop[4])])
  dev.off()
}
## Loop through and make 2D scatter plots for tests
for (i in populations) {
  pdf(file = paste0("Results/ctr 2D scatter/", i ,"_positive_", "_tests.pdf"))
  par(mfrow = c(2,2), pch = '.')
  pop <- populations[populations != i]
  plot(get(i)[,c(i, pop[1])])
  plot(get(i)[,c(i, pop[2])])
  plot(get(i)[,c(i, pop[3])])
  plot(get(i)[,c(i, pop[4])])
  dev.off()
}

## make same scatter plot but put tests and ctr together and with different colours
pdf(file = paste0("Results/ctr 2D scatter/","scatter plots", ".pdf"))
for (i in populations) {
  ctr <- get(paste0(i,"_ctr"))
  ctr$colour <- "#ff000050"
  ctr$Group <- "control"
  test <- get(i)
  test$colour <- "#0000ff50"
  test$Group <- "test"
  data <- rbind(ctr, test)
  #pdf(file = paste0("Results/ctr 2D scatter/", i ,"combined", ".pdf"))
  par(mfrow = c(2,2), pch = 1, cex = 1, xpd = FALSE, bty = "l")
  pop <- populations[populations != i]
  plot(data[,c(i, pop[1])], col = data[,"colour"], main = i)
  plot(data[,c(i, pop[2])], col = data[,"colour"])
  plot(data[,c(i, pop[3])], col = data[,"colour"])
  plot(data[,c(i, pop[4])], col = data[,"colour"])
  legend("topright", c("test", "control"), pch = 1, col = c("#0000ffaa","#ff0000aa"))
  #dev.off()
}
dev.off()
rm(ctr,test,data,pop,i)
} 

## graph intensity distribution of all channels in all populations
## couldnt figure out how to save them, will save manually
for (i in populations) {
  x <- melt(get(i))
  #pdf(file = paste0("Results/", i, "intensity_all_channels_box.pdf"))
  assign(paste0(i,"_boxplot"), ggplot(data = x, aes(x = variable,y = value, fill=variable), legend = FALSE) + geom_boxplot(show.legend = FALSE))
  #dev.off()
}
rm(x)

## Saving boxplots manually
pdf(file = paste0("Results/", "mCherry", "intensity_all_channels_box.pdf" ))
mCherry_boxplot
dev.off()


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

for (i in populations) {
  ## population to be compensated
  target <- paste0(i, "_boxplot")
  pos_median <- median( get(target)$data$value[ get(target)$data$variable == i ] )
  print(pos_median)
  ## go through each channel and calculate k
  for (j in populations) {
    if (j != i) {
      neg_median <- median( get(target)$data$value[ get(target)$data$variable == j ] )
      k <- (neg_median-background[j])/pos_median
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


## Make scatter plots showing each ctr against each others.
combinations <- combn(populations, 2)
colours <- c("#33ccff70", "#00ff0050", "#ffff0050", "#ff990050", "#ff000050")
names(colours) <- populations
pdf(file = "Results/ctr 2D scatter/ctr_comparisons.pdf")
par(pch = 16, cex = 1, xpd = FALSE, bty = "l")
for (i in 1:ncol(combinations)) {
  
  pop_1 <- get(paste0(combinations[1,i], "_ctr"))
  pop_1$colour <- colours[names(colours) == combinations[1,i]]
  pop_2 <- get(paste0(combinations[2,i], "_ctr"))
  pop_2$colour <- colours[names(colours) == combinations[2,i]]
  
  data <- rbind(pop_1, pop_2)
  plot(data[,c(combinations[1,i], combinations[2,i])], col = data$colour, xlim = c(0, 4100), ylim = c(0, 4100))
  legend("topright", c(combinations[1,i], combinations[2,i]), pch = 16, 
         col = c(colours[names(colours) == combinations[1,i]],
                 colours[names(colours) == combinations[2,i]]))
  
}
dev.off()
rm(pop_1, pop_2, data, colours, combinations)

## Do the same but on tests
combinations <- combn(populations, 2)
colours <- c("#33ccff70", "#00ff0050", "#ffff0050", "#ff990050", "#ff000050")
names(colours) <- populations
pdf(file = "Results/ctr 2D scatter/test_comparisons.pdf")
par(pch = 16, cex = 1, xpd = FALSE, bty = "l")
for (i in 1:ncol(combinations)) {
  
  pop_1 <- get(paste0(combinations[1,i]))
  pop_1$colour <- colours[names(colours) == combinations[1,i]]
  pop_2 <- get(paste0(combinations[2,i]))
  pop_2$colour <- colours[names(colours) == combinations[2,i]]
  
  data <- rbind(pop_1, pop_2)
  plot(data[,c(combinations[1,i], combinations[2,i])], col = data$colour, xlim = c(0, 4100), ylim = c(0, 4100))
  legend("topright", c(combinations[1,i], combinations[2,i]), pch = 16, 
         col = c(colours[names(colours) == combinations[1,i]],
                 colours[names(colours) == combinations[2,i]]))
  
}
dev.off()
rm(pop_1, pop_2, data, colours, combinations)

## Make scatter plot of all ctr all shown together and in log scale.
## set y and x max to max in the channel
combinations <- combn(populations, 2)
colours <- c("#33ccff70", "#00ff0050", "#ffff0050", "#ff990050", "#ff000050")
names(colours) <- populations
pdf(file = "Results/ctr 2D scatter/ctr_all_together_comparisons.pdf")
par(pch = 16, cex = 1, xpd = FALSE, bty = "l")
for (i in 1:ncol(combinations)) {
  
  pop_1 <- get(paste0(combinations[1,i], "_ctr"))
  pop_1$colour <- colours[names(colours) == combinations[1,i]]
  pop_1$pop <- combinations[1,i]
  pop_2 <- get(paste0(combinations[2,i], "_ctr"))
  pop_2$colour <- colours[names(colours) == combinations[2,i]]
  pop_2$pop <- combinations[2,i]
  
  neg <- populations[populations != combinations[1,i] & populations != combinations[2,i]]
  pop_3 <- get(paste0(neg[1], "_ctr"))
  pop_3$colour <- colours[names(colours) == neg[1]]
  pop_3$pop <- neg[1]
  pop_4 <- get(paste0(neg[2], "_ctr"))
  pop_4$colour <- colours[names(colours) == neg[2]]
  pop_4$pop <- neg[2]
  pop_5 <- get(paste0(neg[3], "_ctr"))
  pop_5$colour <- colours[names(colours) == neg[3]]
  pop_5$pop <- neg[3]
  
  data <- rbind(pop_1, pop_2, pop_3, pop_4, pop_5)
  data[,1:5] <- data[,1:5] + 1
  data[,1:5] <- log10(data[,1:5])
  plot(data[,c(combinations[1,i], combinations[2,i])], col = data$colour, 
       xlim = c(0
                , data[which.max(data[,combinations[1,i]]), combinations[1,i]]), 
       ylim = c(0
                , data[which.max(data[,combinations[2,i]]), combinations[2,i]])
       
       )
  plot.new()
  legend("topright", c(combinations[1,i], combinations[2,i], neg[1], neg[2], neg[3]), pch = 16, 
         col = c(colours[names(colours) == combinations[1,i]],
                 colours[names(colours) == combinations[2,i]],
                 colours[names(colours) == neg[1]],
                 colours[names(colours) == neg[2]],
                 colours[names(colours) == neg[3]]))
  
}
dev.off()
rm(pop_1, pop_2, pop_3, pop_4, pop_5, data, colours, combinations)

## Graph put tests on to the last graph
combinations <- combn(populations, 2)
colours <- c("#33ccff40", "#00ff0020", "#ffff0020", "#ff990020", "#ff000020")
names(colours) <- populations
pdf(file = "Results/ctr 2D scatter/ctr_test_together_comparisons.pdf")
par( xpd = FALSE, bty = "l")
for (i in 1:ncol(combinations)) {
  
  colours <- c("#33ccff50", "#00ff0050", "#e6e60050", "#ff990050", "#ff000050")
  names(colours) <- populations
  pop_1 <- get(paste0(combinations[1,i], "_ctr"))
  pop_1$colour <- colours[names(colours) == combinations[1,i]]
  pop_1$pop <- combinations[1,i]
  pop_2 <- get(paste0(combinations[2,i], "_ctr"))
  pop_2$colour <- colours[names(colours) == combinations[2,i]]
  pop_2$pop <- combinations[2,i]
  
  neg <- populations[populations != combinations[1,i] & populations != combinations[2,i]]
  pop_3 <- get(paste0(neg[1], "_ctr"))
  pop_3$colour <- colours[names(colours) == neg[1]]
  pop_3$pop <- neg[1]
  pop_4 <- get(paste0(neg[2], "_ctr"))
  pop_4$colour <- colours[names(colours) == neg[2]]
  pop_4$pop <- neg[2]
  pop_5 <- get(paste0(neg[3], "_ctr"))
  pop_5$colour <- colours[names(colours) == neg[3]]
  pop_5$pop <- neg[3]
  ctr <- rbind(pop_1, pop_2, pop_3, pop_4, pop_5)
  ctr$pch <- 18
  
  colours <- c("#33ccffaa", "#00ff00aa", "#e6e600aa", "#ff9900aa", "#ff0000aa")
  names(colours) <- populations
  ## read in test populations and give it pch value and colour values
  for (j in 1:length(populations)) {
    assign(paste0("test_", j), get(populations[j]))
  }
  test_1$colour <- colours[names(colours) == populations[1]]
  test_1$pop <- populations[1]
  test_2$colour <- colours[names(colours) == populations[2]]
  test_2$pop <- populations[2]
  test_3$colour <- colours[names(colours) == populations[3]]
  test_3$pop <- populations[3]
  test_4$colour <- colours[names(colours) == populations[4]]
  test_4$pop <- populations[4]
  test_5$colour <- colours[names(colours) == populations[5]]
  test_5$pop <- populations[5]
  test <- rbind(test_1, test_2, test_3, test_4, test_5)
  test$pch <- 21
  
  data <- rbind(ctr, test)
  data[,1:5] <- data[,1:5] + 1
  data[,1:5] <- log10(data[,1:5])
  plot(data[,c(combinations[1,i], combinations[2,i])], col = data$colour, 
       xlim = c(0
                , data[which.max(data[,combinations[1,i]]), combinations[1,i]]), 
       ylim = c(0
                , data[which.max(data[,combinations[2,i]]), combinations[2,i]]),
       pch = data$pch
  )
  plot.new()
  legend("topright", c("CFP ctr", "GFP ctr", "YFP ctr", "mKO2 ctr", "mCherry ctr"), pch = 18, 
         col = c(colours))
  legend("right", c("CFP test", "GFP test", "YFP test", "mKO2 test", "mCherry test"), pch = 21, 
         col = colours
         )
}
dev.off()
rm(pop_1, pop_2, pop_3, pop_4, pop_5, data, colours, combinations)






plot <- ggplot(data = mKO2, aes(x = mKO2, y = YFP)) + geom_point(shape = 1)
plot

x <- melt(YFP.comp[,c(1,3)])
ggplot(data = x,aes(x = value, fill = variable)) + geom_density(alpha=0.25)
ggplot(data = x,aes(x = variable,y = value, fill=variable)) + geom_boxplot(show.legend = FALSE)


