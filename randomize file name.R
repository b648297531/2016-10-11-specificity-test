
## get file name and assign random number to them
sample_list <- list.files(path = "Test images copy/", pattern = ".nd2")

names(sample_list) <- sample(c(1:25), size = 25)

## change file names according to its name in the list

for (i in 1:25) {
  file.rename(from = paste0("Test images copy/", sample_list[i]), to = paste0("Test images copy/", names(sample_list[i]), ".nd2"))
}
write.table(sample_list, file = "table.txt")


