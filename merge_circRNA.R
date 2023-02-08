library(tidyverse)
all_files <- dir("/home/suqiang/90-samples/bam_reads/circRNA_finder/filteredJun/csv_data")
file_names <- grep(all_files,pattern = "*99272.csv",value = TRUE)
#file_names
Data_file <- map(file_names,read.delim,header = FALSE, stringsAsFactors = FALSE, check.names = FALSE, row.names = NULL, sep = ' ')
#Data_file
Merge_All_Samples <- Data_file %>% reduce(inner_join, by = "V1")
write.csv(Merge_All_Samples, "./Merge_All_N_Samples.csv", row.names = F)





