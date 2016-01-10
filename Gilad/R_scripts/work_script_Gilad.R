# R commands to reanalyze the data from Lin et al (2014) starting from the raw fragment counts

#0. Set working directory, change as needed
setwd("/Users/mikhail/Documents/Work/GenomeRunner/bodymap/Gilad/R_input_files")

#1. load data
datasets = as.data.frame(scan("Stanford_datasets.txt",list(setname="",seqBatch="",species="",tissue=""),sep="\t"))
rawCounts <- as.matrix(read.table('Stanford_datasets_rawCountsMat.txt',header=FALSE,sep='\t'))
colnames(rawCounts) <- datasets$setname
geneDetails <- as.data.frame(scan("ortholog_GC_table.txt",skip=1,list(mouse_name="",mouse_GC = 0.0,human_name = "",human_GC=0.0)))
rownames(rawCounts) <- geneDetails$human_name 
rownames(datasets) <- datasets$setname
rownames(geneDetails) <- geneDetails$human_name 

#2. filter out lower 30% + mitochondrial genes
rowSums = apply(rawCounts,1,function(x) sum(x))
quantile(rowSums,probs = 0.3) # result is 2947.9 
filter <- apply(rawCounts,1,function(x) sum(x)>2947.9 )
mt <- grep("mt-",geneDetails$mouse_name)
filteredNames <- setdiff(rownames(rawCounts[filter,]),rownames(rawCounts[mt,])) 
filteredRawCounts <- rawCounts[filteredNames,]

#3. Use EDASeq to normalize data within lanes, accounting for GC content
library(EDASeq)
GCnormCounts <- filteredRawCounts
GCnormCounts[,1:13] <- withinLaneNormalization(filteredRawCounts[,1:13],geneDetails[filteredNames,"human_GC"],which="loess",round=TRUE)
GCnormCounts[,14:26] <- withinLaneNormalization(filteredRawCounts[,14:26],geneDetails[filteredNames,"mouse_GC"],which="loess",round=TRUE)

#4. depth normalize,using TMM scaling factors - divide by sum, then multiply by mean of sums
library(edgeR)
origColSums <- apply(rawCounts,2,function(x) sum(x))
normFactors <- calcNormFactors(GCnormCounts,method='TMM')
normalizedColSums <- origColSums
i <- 1
while (i<length(normalizedColSums)){
  normalizedColSums[i] <- origColSums[i]* normFactors[i]
  i <- i+1
}
meanDepth <- mean(normalizedColSums)
filteredDepthNormCounts <- GCnormCounts
i <- 1
while (i<ncol(filteredDepthNormCounts)){
  filteredDepthNormCounts[,i] <- (GCnormCounts[,i]/normalizedColSums[i])*meanDepth
  i <- i+1
}

#5. log transformation
logTransformedDepthNormCounts <- log2(filteredDepthNormCounts+1)

#6. use combat on log2 values to remove batch effects 
library(sva)
meta <- data.frame(seqBatch = datasets$seqBatch,tissue=datasets$tissue,species=datasets$species)
design <- model.matrix(~1,data=meta)
combat <- ComBat(dat= logTransformedDepthNormCounts,batch=meta$seqBatch,mod=design,par.prior=TRUE)

#7. Save the depth-normalized, GC- and batch adjusted log2-transformed counts data
write.table(combat, "Stanford_dataset_normalized_counts_full.txt", sep = "\t", quote = F, col.names = NA)
