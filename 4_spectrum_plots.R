## plot spectra for ctrs

## read files


populations <- c("GFP", "YFP", "mKO2", "mCherry")

## read data
for (i in populations) {
  #assign(i, read.table(file = paste0("Test images (individual)/Measurements/",i,".txt"))) ## test images
  assign(paste0(i,"_ctr"), read.table(file = paste0("unmix CTR/Measurements/",i,".txt"))) ## ctr images
}

for (i in populations) {
  
  file <- get(paste0(i,"_ctr"))
  rownames(file) <- 1:nrow(file)
  for (j in 1:ncol(file)) {
    colnames(file)[j] <- 434+10*(j-1)
  }
  colnames(file)[ncol(file)] <- "TD"
  assign(paste0(i,"_ctr"), file)
}
rm(file, i, j)



## plot scatter plot for all
### y-axis: intensity
### x axis: wavelength



for (i in populations) {
  data <- matrix(nrow = ncol(get(paste0(i, "_ctr")))*nrow(get(paste0(i, "_ctr"))), ncol = 2)
  a <- 1
  for (j in 1:(ncol(get(paste0(i, "_ctr")))-1)) {
    data[a:(nrow(get(paste0(i, "_ctr")))*j),1] <- colnames(get(paste0(i, "_ctr")))[j]
    data[a:(nrow(get(paste0(i, "_ctr")))*j),2] <- get(paste0(i, "_ctr"))[,j]
    a <- nrow(get(paste0(i, "_ctr")))*j+1
    print(a)
  }
  data[,2] <- as.numeric(data[,2])/as.numeric(data[which.max(data[,2]),2])
  
  
  pdf(file = paste0("Results/Spectra/", i, "_ctr.pdf"))
  plot(data, pch =16, bty = "l", xlab = "Wavelenth (nm)", ylab = "Fluorescence intensity")
  dev.off()
}


## calculate avarage in each channel and plot spectra

pdf(file = paste0("Results/Spectra/all_ctr_mean_median.pdf"))
for (i in populations) {
  
  # wavelength, average intensity, 1 sd up, 1 sd down, median, 1st Qu, 3rd Qu
  data <- matrix(ncol = 7, nrow = ncol(get(paste0(i, "_ctr")))-1 )

  for (j in 1:( ncol(get(paste0(i, "_ctr")))-1 ) ) {
    data[j,1] <- as.numeric(colnames(get(paste0(i, "_ctr")))[j] )
    data[j,2] <- mean( get(paste0(i, "_ctr"))[,j] / max(get(paste0(i, "_ctr"))[,1:32]) )
    data[j,3] <- data[j,2] + sd(get(paste0(i, "_ctr"))[,j] / max(get(paste0(i, "_ctr"))[,1:32] ))
    data[j,4] <- data[j,2] - sd(get(paste0(i, "_ctr"))[,j] / max(get(paste0(i, "_ctr"))[,1:32] ))
    data[j,5] <- median( get(paste0(i, "_ctr"))[,j] / max(get(paste0(i, "_ctr"))[,1:32]) )
    data[j,6] <- summary(get(paste0(i, "_ctr"))[,j] / max(get(paste0(i, "_ctr"))[,1:32]) )[2] 
    data[j,7] <- summary(get(paste0(i, "_ctr"))[,j] / max(get(paste0(i, "_ctr"))[,1:32]) )[5] 
  }
   
  #pdf(file = paste0("Results/Spectra/", i, "_ctr_mean.pdf"))
  plot(data[,1:2], pch =16, bty = "l", xlab = "Wavelenth (nm)", 
       ylab = paste0("Mean ", i, " fluorescence intensities relative to maximum"),
       ylim = c(0, data[which.max(data[,3]),3]), las = 1)
  lines(data[,1:2])
  arrows(data[,1], data[,3], data[,1], data[,4], code = 3, length = 0.05, angle = 90, col = "red")
  #dev.off()
  
  #pdf(file = paste0("Results/Spectra/", i, "_ctr_median.pdf"))
  plot(data[,c(1,5)], pch =16, bty = "l", xlab = "Wavelenth (nm)", 
       ylab = paste0("Median ",i, " fluorescence intensities relative to maximum"),
       ylim = c(0, data[which.max(data[,3]),3]), las = 1)
  lines(data[,c(1,5)])
  arrows(data[,1], data[,6], data[,1], data[,7], code = 3, length = 0.05, angle = 90, col = "red")
  #dev.off()
  
  data <- as.data.frame(data)
  colnames(data) <- c("wavelength", "mean", "sd.up", "sd.down", "median", "1stQu", "3rdQu")
  assign(paste0(i, "_ctr_mean"), data)
  
}
dev.off()

write.table(GFP_ctr_mean, file = "Results/Spectra/GFP_ctr_sum.txt", sep = '\t')
write.table(YFP_ctr_mean, file = "Results/Spectra/YFP_ctr_sum.txt", sep = '\t')
write.table(mKO2_ctr_mean, file = "Results/Spectra/mKO2_ctr_sum.txt", sep = '\t')
write.table(mCherry_ctr_mean, file = "Results/Spectra/mCherry_ctr_sum.txt", sep = '\t')


#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
######### Do same things for the test populations ###############################################################


## read files

populations <- c("GFP", "YFP", "mKO2", "mCherry")

for (i in populations) {
  #assign(i, read.table(file = paste0("Test images (individual)/Measurements/",i,".txt"))) ## test images
  assign(paste0(i), read.table(file = paste0("Test images/Measurements/",i,".txt"))) ## ctr images
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


## calculate avarage in each channel and plot spectra

pdf(file = paste0("Results/Spectra/all_test_mean_median.pdf"))
for (i in populations) {
  
  # wavelength, average intensity, 1 sd up, 1 sd down, median, 1st Qu, 3rd Qu
  data <- matrix(ncol = 7, nrow = ncol(get(i))-1 )
  
  for (j in 1:( ncol(get(i))-1 ) ) {
    data[j,1] <- as.numeric(colnames(get(i))[j] )
    data[j,2] <- mean( get(i)[,j] / max(get(i)[,1:32]) )
    data[j,3] <- data[j,2] + sd(get(i)[,j] / max(get(i)[,1:32] ))
    data[j,4] <- data[j,2] - sd(get(i)[,j] / max(get(i)[,1:32] ))
    data[j,5] <- median( get(i)[,j] / max(get(i)[,1:32]) )
    data[j,6] <- summary(get(i)[,j] / max(get(i)[,1:32]) )[2] 
    data[j,7] <- summary(get(i)[,j] / max(get(i)[,1:32]) )[5] 
  }
  
  #pdf(file = paste0("Results/Spectra/", i, "_ctr_mean.pdf"))
  plot(data[,1:2], pch =16, bty = "l", xlab = "Wavelenth (nm)", 
       ylab = paste0("Mean ", i, " fluorescence intensities relative to maximum"),
       ylim = c(0, data[which.max(data[,3]),3]), las = 1)
  lines(data[,1:2])
  arrows(data[,1], data[,3], data[,1], data[,4], code = 3, length = 0.05, angle = 90, col = "red")
  #dev.off()
  
  #pdf(file = paste0("Results/Spectra/", i, "_ctr_median.pdf"))
  plot(data[,c(1,5)], pch =16, bty = "l", xlab = "Wavelenth (nm)", 
       ylab = paste0("Median ",i, " fluorescence intensities relative to maximum"),
       ylim = c(0, data[which.max(data[,3]),3]), las = 1)
  lines(data[,c(1,5)])
  arrows(data[,1], data[,6], data[,1], data[,7], code = 3, length = 0.05, angle = 90, col = "red")
  #dev.off()
  
  data <- as.data.frame(data)
  colnames(data) <- c("wavelength", "mean", "sd.up", "sd.down", "median", "1stQu", "3rdQu")
  assign(paste0(i, "_test_mean"), data)
  
}
dev.off()



## plot test vs ctr in same graph

pdf(file = paste0("Results/Spectra/all_ctr_test_mean_median.pdf"))
for (i in populations) {
  
  # wavelength, average intensity, 1 sd up, 1 sd down, median, 1st Qu, 3rd Qu, colour
  data_test <- matrix(ncol = 8, nrow = ncol(get(paste0(i, "_ctr")))-1 )
  
  for (j in 1:( ncol(get(i))-1 ) ) {
    data_test[j,1] <- as.numeric(colnames(get(i))[j] )
    data_test[j,2] <- mean( get(i)[,j] )
    data_test[j,3] <- data_test[j,2] + sd(get(i)[,j])
    data_test[j,4] <- data_test[j,2] - sd(get(i)[,j])
    data_test[j,5] <- median( get(i)[,j] )
    data_test[j,6] <- summary(get(i)[,j])[2] 
    data_test[j,7] <- summary(get(i)[,j])[5] 
  }
  data_test[,8] <- "blue"
  
  data_ctr <- matrix(ncol = 8, nrow = ncol(get(paste0(i, "_ctr")))-1 )
  
  for (j in 1:( ncol(get(paste0(i, "_ctr")))-1 ) ) {
    data_ctr[j,1] <- as.numeric(colnames(get(paste0(i, "_ctr")))[j] )
    data_ctr[j,2] <- mean( get(paste0(i, "_ctr"))[,j])
    data_ctr[j,3] <- data_ctr[j,2] + sd(get(paste0(i, "_ctr"))[,j])
    data_ctr[j,4] <- data_ctr[j,2] - sd(get(paste0(i, "_ctr"))[,j])
    data_ctr[j,5] <- median( get(paste0(i, "_ctr"))[,j] )
    data_ctr[j,6] <- summary(get(paste0(i, "_ctr"))[,j])[2] 
    data_ctr[j,7] <- summary(get(paste0(i, "_ctr"))[,j])[5]
  }
  data_ctr[,8] <- "red"
  
  data <- rbind(data_ctr, data_test)
  #pdf(file = paste0("Results/Spectra/", i, "_ctr_mean.pdf"))
  plot(data[,1:2] , pch =16, bty = "l", xlab = "Wavelenth (nm)", 
       ylab = paste0("Mean ", i, " fluorescence intensities"),
       ylim = c(0, as.numeric(data[which.max(data[,3]),3])), las = 1, 
       col = data[,8])
  lines(data[data[,8] == "red",1:2])
  lines(data[data[,8] == "blue",1:2])
  arrows(as.numeric(data[,1]), 
         as.numeric(data[,3]), 
         as.numeric(data[,1]), 
         as.numeric(data[,4]), code = 3, length = 0.05, angle = 90, col = data[,8])
  legend("topright", c("Control population", "Test population"),
         pch = 16, 
         col = c("red", "blue")
         ,cex = 1)
  #dev.off()
  
  #pdf(file = paste0("Results/Spectra/", i, "_ctr_median.pdf"))
  plot(data[,c(1,5)], pch =16, bty = "l", xlab = "Wavelenth (nm)", 
       ylab = paste0("Median ", i, " fluorescence intensities"),
       ylim = c(0, as.numeric(data[which.max(data[,7]),7])), las = 1, 
       col = data[,8])
  lines(data[data[,8] == "red",c(1,5)])
  lines(data[data[,8] == "blue",c(1,5)])
  arrows(as.numeric(data[,1]), 
         as.numeric(data[,6]), 
         as.numeric(data[,1]), 
         as.numeric(data[,7]), code = 3, length = 0.05, angle = 90, col = data[,8])  #dev.off()
  legend("topright", c("Control population", "Test population"),
         pch = 16, 
         col = c("red", "blue")
         ,cex = 1)
  
}
dev.off()




#######################################################################################################################
#######################################################################################################################
## plot all ctr pop in one graph

colours <- c( "#00ff00a0", "#ffff00ff", "#ff9900a0", "#ff0000a0")
names(colours) <- populations
points <- c(21,22,23,24)
names(points) <- populations

data <- matrix(ncol = 10, nrow = 1)
for(i in populations) {
  
  data_ctr <- matrix(ncol = 10, nrow = ncol(get(paste0(i, "_ctr")))-1 )
  for (j in 1:( ncol(get(paste0(i, "_ctr")))-1 ) ) {
    data_ctr[j,1] <- as.numeric(colnames(get(paste0(i, "_ctr")))[j] )
    data_ctr[j,2] <- mean( get(paste0(i, "_ctr"))[,j])
    data_ctr[j,3] <- data_ctr[j,2] + sd(get(paste0(i, "_ctr"))[,j])
    data_ctr[j,4] <- data_ctr[j,2] - sd(get(paste0(i, "_ctr"))[,j])
    data_ctr[j,5] <- median( get(paste0(i, "_ctr"))[,j] )
    data_ctr[j,6] <- summary(get(paste0(i, "_ctr"))[,j])[2] 
    data_ctr[j,7] <- summary(get(paste0(i, "_ctr"))[,j])[5]
  }
  data_ctr[,9] <- i
  data_ctr[,8] <- colours[names(colours) == i]
  data_ctr[,10] <- points[names(points) == i]
  
  data <- rbind(data_ctr, data)
}
data <- data[1:(nrow(data)-1),]

pdf(file = paste0("Results/Spectra/all_ctr_combined_mean_median.pdf"))
plot(data[,c(1,5)], pch = as.numeric(data[,10]), bty = "l", xlab = "Wavelenth (nm)", 
     ylab = paste0("Median ", i, " fluorescence intensities"),
     ylim = c(0, as.numeric(data[which.max(data[,7]),7])), las = 1, 
     col = data[,8])
for (i in populations) {
  lines(data[data[,9] == i,c(1,5)], col = colours[names(colours) == i])
}

arrows(as.numeric(data[,1]), 
       as.numeric(data[,6]), 
       as.numeric(data[,1]), 
       as.numeric(data[,7]), code = 3, length = 0.05, angle = 90, col = data[,8])  #dev.off()

legend("topright", populations,
       pch = points, 
       col = colours,
       cex = 1)
dev.off()


####################################################################################################################################
####################################################################################################################################
###### plot all tests in one graph
colours <- c( "#00ff00a0", "#ffff00ff", "#ff9900a0", "#ff0000a0")
names(colours) <- populations
points <- c(21,22,23,24)
names(points) <- populations

data <- matrix(ncol = 10, nrow = 1)
for(i in populations) {
  
  data_ctr <- matrix(ncol = 10, nrow = ncol(get(i))-1 )
  for (j in 1:( ncol(get(i))-1 ) ) {
    data_ctr[j,1] <- as.numeric(colnames(get(i))[j] )
    data_ctr[j,2] <- mean( get(i)[,j])
    data_ctr[j,3] <- data_ctr[j,2] + sd(get(i)[,j])
    data_ctr[j,4] <- data_ctr[j,2] - sd(get(i)[,j])
    data_ctr[j,5] <- median( get(i)[,j] )
    data_ctr[j,6] <- summary(get(i)[,j])[2] 
    data_ctr[j,7] <- summary(get(i)[,j])[5]
  }
  data_ctr[,9] <- i
  data_ctr[,8] <- colours[names(colours) == i]
  data_ctr[,10] <- points[names(points) == i]
  
  data <- rbind(data_ctr, data)
}
data <- data[1:(nrow(data)-1),]

pdf(file = paste0("Results/Spectra/all_test_combined_mean_median.pdf"))
plot(data[,c(1,5)], pch = as.numeric(data[,10]), bty = "l", xlab = "Wavelenth (nm)", 
     ylab = paste0("Median ", i, " fluorescence intensities"),
     ylim = c(0, as.numeric(data[which.max(data[,7]),7])), las = 1, 
     col = data[,8])
for (i in populations) {
  lines(data[data[,9] == i,c(1,5)], col = colours[names(colours) == i])
}

arrows(as.numeric(data[,1]), 
       as.numeric(data[,6]), 
       as.numeric(data[,1]), 
       as.numeric(data[,7]), code = 3, length = 0.05, angle = 90, col = data[,8])  #dev.off()

legend("topright", populations,
       pch = points, 
       col = colours,
       cex = 1)
dev.off()












































































