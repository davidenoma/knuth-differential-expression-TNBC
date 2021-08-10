# #Install packages
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.13")
# BiocManager::install("DESeq2")
# BiocManager::install("clusterProfiler")
# BiocManager::install("org.Hs.eg.db")
# BiocManager::install("apeglm")
# BiocManager::install("EnsDb.Hsapiens.v86")
# BiocManager::install("pheatmap")
# BiocManager::install("biomaRt")
# BiocManager::install("RColorBrewer")
# BiocManager::install("GEOquery")
#load packages
library(DESeq2)

library(clusterProfiler)
library(org.Hs.eg.db)
library(apeglm)
library(EnsDb.Hsapiens.v86)
library(pheatmap)
library(pheatmap)
library(RColorBrewer)
library(biomaRt)
library(GEOquery)
library(clusterProfiler)

gse142731 <- GEOquery::getGEO("GSE142731")
print(gse142731)
colnames(pData(gse142731[[1]]))[1:10]
pheno_data=pData(gse142731[[1]])


#Quick exploratory analysis
print(table(pheno_data$`race:ch1`))
print(table(pheno_data$characteristics_ch1))
print(table(pheno_data$`tissue:ch1`))

#declare a dataframe
labels_gse142731=data.frame(V1=pData(gse142731[[1]])[,"race:ch1"])
#convert into a factor variable
labels_gse142731$V1=factor(labels_gse142731$V1)
levels(labels_gse142731$V1)=c('African American','Caucasian')

GSE142731<-as.matrix(read.table("C:/Users/HP/Downloads/GSE142731_Matrix-file.txt"))
dim(GSE142731)
deseq2_142731 <- DESeqDataSetFromMatrix(countData = round(GSE142731[,2:ncol(GSE142731)]),colData = labels_gse142731,design = ~V1)

#creating a copy of the DESeq2 object
dds_142731=deseq2_142731 
#making a copy of this DESeq2 object to get a list of DEGs


#Removing rows with zero values form the dataset
keep_genes <- rowSums(counts(deseq2_142731)) > 0
deseq2_142731 <- deseq2_142731[ keep_genes, ]
dim(deseq2_142731)



#size factor estimation
deseq2_142731 <- estimateSizeFactors(deseq2_142731)
#difference between non-normalized and log-normalized read counts
#Log transformed boxlpots
boxplot(log2(counts(deseq2_142731)+1), notch=TRUE,
        main = "Non-normalized read counts",
        ylab="log2(read counts)", cex = .6,xaxt="n")
## bp of size-factor normalized values
boxplot(log2(counts(deseq2_142731, normalize= TRUE) +1), notch=TRUE,
        main = "Size-factor-normalized read counts",
        ylab="log2(read counts)", cex = .6,xaxt="n")


#VST
vsd_142731 <- vst(deseq2_142731, blind = FALSE)
head(assay(vsd_142731), 3)
#rlog
rld <- rlog(deseq2_142731, blind = FALSE)
rld_142731=rld
head(assay(rld_142731), 3)
#plot the boxplots of the expression measurements for all 3 methods
par(mfrow=c(1,3))
plot(assay(vsd_142731),main="VST")
plot(assay(rld_142731),main="rlog")
plot(log2(counts(deseq2_142731, normalized=TRUE)[, 1:2]+1),main="log2")


#CLuster samples across all genes
sampleDists_142731 <- dist(t(assay(vsd_142731)))
sampleDists_142731
sampleDistMatrix_142731 <- as.matrix( sampleDists_142731 )
rownames(sampleDistMatrix_142731) <- as.character(vsd_142731$V1)
colnames(sampleDistMatrix_142731) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix_142731,
         clustering_distance_rows = sampleDists_142731,
         clustering_distance_cols = sampleDists_142731,
         col = colors)
#PCA plot
plotPCA(vsd_142731, intgroup = c("V1"))


#Differential expression analysis
dds_142731$V1 <- relevel(dds_142731$V1, ref = "African American")
dds_142731 <- DESeq(dds_142731)

resultsNames(dds_142731)

res_142731=results(dds_142731,contrast=c("V1","Caucasian","African American"))
res_142731 = na.omit(res_142731) #removing NA
summary(res_142731)
plotMA(res_142731)

#FDR cutoff 0.05
resfiltered_142731 = res_142731[res_142731$padj <= 0.05,]
#order by adjusted p-value first
res.filtered.ordered_142731 = 
  resfiltered_142731[order(resfiltered_142731$padj),] 
#convert to a dataframe
res.filtered.ordered_142731=as.data.frame(res.filtered.ordered_142731)
#Seperate up regulated genes
Sign_genes_up_142731 = subset(res.filtered.ordered_142731, log2FoldChange >0.5)
#Seperate up regulated genes
Sign_genes_down_142731 = subset(res.filtered.ordered_142731, log2FoldChange < -0.5)
#make a consolidated list with up-regulated genes in the begininig
genes_142731=c(rownames(Sign_genes_down_142731),rownames(Sign_genes_up_142731))

top10_down<-head(Sign_genes_down_142731[order(Sign_genes_down_142731$padj),], 10)
top10_up<-head(Sign_genes_down_142731[order(Sign_genes_down_142731$padj),], 10)



#Gene Ontology Analysis
all_genes_human = unique(rownames(dds_142731))
select(org.Hs.eg.db,rownames(Sign_genes_up_142731),c("SYMBOL","GENENAME","ENTREZID"),"ENSEMBL")
select(org.Hs.eg.db,rownames(Sign_genes_down_142731),c("SYMBOL","GENENAME","ENTREZID"),"ENSEMBL")
select(org.Hs.eg.db,rownames(top10_down),c("SYMBOL","GENENAME","ENTREZID"),"ENSEMBL")
select(org.Hs.eg.db,rownames(top10_up),c("SYMBOL","GENENAME","ENTREZID"),"ENSEMBL")

