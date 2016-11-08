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

## For measurements on spectral imaging of ctr

populations <- c("GFP", "YFP", "mKO2", "mCherry") #* CFP is dropped

### read all the files. 
for (i in populations) {
    file_name <- paste0(i, ".xls")
    file <- t(read.table(file = paste0("unmix CTR/Measurements/", file_name), 
                         sep = '\t', header = TRUE, row.names = 1))
    rownames(file) <- 1:nrow(file)
    for (j in 1:ncol(file)) {
      colnames(file)[j] <- 434+10*(j-1)
    }
    colnames(file)[ncol(file)] <- "TD"
    assign(paste0(i), value = as.data.frame(file))
}
rm(i,j,file, file_name)

for (i in populations) {
  write.table(get(i), file = paste0("unmix CTR/Measurements/", i, ".txt"), sep = '\t')
}


## For measurements on spectal imaging os test
populations <- c("CFP", "GFP", "YFP", "mKO2", "mCherry")

### read all the files. 
for (i in populations) {
  for (j in 1:5){
    file_name <- paste0(i, j, ".xls")
    file <- t(read.table(file = paste0("Test images/Measurements/", file_name), 
                         sep = '\t', header = TRUE, row.names = 1))
    rownames(file) <- 1:nrow(file)
    for (k in 1:ncol(file)) {
      colnames(file)[k] <- 434+10*(k-1)
    }
    colnames(file)[ncol(file)] <- "TD"
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



## Do not need to run those anymore ^


