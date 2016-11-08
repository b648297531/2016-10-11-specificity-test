
library(ggplot2)
library(reshape2)

populations <- c("CFP", "GFP", "YFP", "mKO2", "mCherry")

## read data
for (i in populations) {
  assign(i, read.table(file = paste0("Test images (individual)/Measurements/",i,".txt"))) ## test images
  assign(paste0(i,"_ctr"), read.table(file = paste0("unmixed/Measurements/",i,".txt"))) ## ctr images
}

#########################################################################################################
#########################################################################################################


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

################################################################################################
################################################################################################

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














