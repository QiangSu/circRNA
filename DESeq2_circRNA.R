library(DESeq2)
library(dplyr)
#loading counts matrix
count <- read.csv("J:/shenzhen Uui/project/circular RNA/code/Merge_all_samples-STAR-circRNA_finder_count-N-T.csv", header=T,row.names=1,check.names = F)
#print(count)

## filtering the counts
count <- count[rowMeans(count)>1,]

##loading sample information
data <- read.table("J:/shenzhen Uui/project/circular RNA/code/sample_list-N-T.txt",header = T,row.names = 1)

#setting factor
data[,1] <- as.factor(data$Type)

all(rownames(data) %in% colnames(count))

#all(rownames(data) == colnames(count))

dds <-  DESeqDataSetFromMatrix(countData = count,colData = data,design = ~ Type)
dim(dds)

#过滤
dds <- dds[rowSums(counts(dds)) > 1,]
nrow(dds) 

## processing
dep <- DESeq(dds)
res <- results(dep)
diff = res
diff <- na.omit(diff)  
dim(diff)
#diff
write.csv(diff,"J:/shenzhen Uui/project/circular RNA/code/Merge_all_samples-STAR-circRNA_finder_count-all_diff_N-T.csv")

#foldChange = 1
#padj = 0.05

#diffsig <- diff[(diff$pvalue < padj & abs(diff$log2FoldChange) > foldChange),]
#dim(diffsig)

#write.csv(diffsig, "J:/shenzhen Uui/project/circular RNA/code/m-merge_all_sample_BWA-CIRI2_count-0131_all_diffsig-N-T.csv")
