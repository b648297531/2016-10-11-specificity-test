## Rewrite 5_criteria_least_distance.R because that was ugly. 
## Make a heat map showing distance and cluster of different populations identified with 5_criteria_least_distance.R


### Import data sets. 
#### Test pop
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

populations <- c("GFP", "YFP", "mKO2", "mCherry")

#### ctr pop
for (i in populations) {
  assign(paste0(i,"_ctr"), read.table(file = paste0("unmixing ctr/Measurements/",i,".txt"))) ## ctr images
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


### Calculating median of ctr populations 
### Save in a table 
###     - Column : median in each wavelength
###     - Row : population identifier 

ctr_medians <- matrix(nrow = 4, ncol = 33)
ctr_medians <- as.data.frame(ctr_medians)
colnames(ctr_medians) <- c(seq(from = 434, to = 744, by = 10), "TD")
rownames(ctr_medians) <- populations

for (i in populations) {
  
  for (j in 1:33) {
    ctr_medians[i,j] <- median(as.numeric( get( paste0(i,"_ctr") )[,j] ))
  }
  
}
rm(i, j)
write.table(ctr_medians, file = "unmixing ctr/Measurements/ctr_medians.txt", sep = '\t')


### Calculate distance of each test cells to ctr median
### Merge all test populations to one data frame with identifier 

GFP$ID <- "GFP"
YFP$ID <- "YFP"
mCherry$ID <- "mCherry"
mKO2$ID <- "mKO2"
Test_pop <- rbind(GFP, YFP, mCherry, mKO2)
rm(GFP, YFP, mCherry, mKO2)
write.table(Test_pop, file = "Test images/Measurements/pooled.txt", sep = '\t')

#### Make table to save distance data
dis_table <- matrix(nrow = nrow(Test_pop), ncol = 5)
dis_table <- as.data.frame(dis_table)
colnames(dis_table) <- c("ID", populations)

for (i in 1:nrow(Test_pop)) {
  dis_table$ID[i] <- Test_pop$ID[i]
  
  for (j in populations) {
    dis_table[i, j] <- dist(rbind(Test_pop[i,1:32]/Test_pop[i,which.max(Test_pop[i,1:32])], 
                                  ctr_medians[j,1:32]/ctr_medians[j,which.max(ctr_medians[j,1:32])]), 
                            method = "euclidean")
  }
  
}
rm(i,j)
write.table(dis_table, file = "Test images/Measurements/euclidean_distance_to_ctr_median(normalized_to_max).txt", sep = '\t')


### Count cells
### Make table. row : cell true ID
###                   cell ID based on least distance
Counts <- data.frame(c(0,0,0,0),
                     c(0,0,0,0),
                     c(0,0,0,0),
                     c(0,0,0,0))
colnames(Counts) <- populations
rownames(Counts) <- populations

for (i in 1:nrow(dis_table)) {
  pos <- names(which.min(dis_table[i, 2:5] ) )
  Counts[dis_table$ID[i], pos] <- Counts[dis_table$ID[i], pos] + 1
}
write.table(Counts, file = "Results/count_least_euclidean_distance_mean.txt", sep = '\t')



####################### Actual heatmapping

#### Read data because I deleted it before
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

#### Add fluorescence ID and test/ctr grouping
GFP$ID <- "GFP"
YFP$ID <- "YFP"
mCherry$ID <- "mCherry"
mKO2$ID <- "mKO2"

GFP_ctr$ID <- "GFP"
YFP_ctr$ID <- "YFP"
mCherry_ctr$ID <- "mCherry"
mKO2_ctr$ID <- "mKO2"

GFP$pop <- "Test"
YFP$pop <- "Test"
mCherry$pop <- "Test"
mKO2$pop <- "Test"

GFP_ctr$pop <- "CTR"
YFP_ctr$pop <- "CTR"
mCherry_ctr$pop <- "CTR"
mKO2_ctr$pop <- "CTR"

#### rbind everything so that i can make the heat map
#### Normalize heatmap data the same way as criteria
heatmap_data <- rbind(GFP, YFP, mCherry, mKO2, GFP_ctr, YFP_ctr, mCherry_ctr, mKO2_ctr)
for (i in 1:nrow(heatmap_data)) {
  
  heatmap_data[i,1:32] <- heatmap_data[i,1:32]/heatmap_data[i,which.max(heatmap_data[i,1:32])]
  
}


library(pheatmap)
annotation_col <- data.frame(Fluorescence=heatmap_data$ID, Population=heatmap_data$pop)


pheatmap(t(heatmap_data[,1:32]), annotation_col = annotation_col, treeheight_col = 100, cluster_rows = FALSE, scale = "column",clustering_method = "median")















































