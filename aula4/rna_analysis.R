library(dplyr)
library(tibble)
##########################################
## Build count matrix with HTSeq counts ##
##########################################

count.files<-c("La_A_1.counts","La_A_2.counts","La_A_3.counts","La_B_1.counts","La_B_2.counts","La_B_3.counts")

# Build the count matrix
count_matrix<-lapply(1:length(count.files), function(x){
  count.vec<-read.delim(file=paste(count.files[x]), header=FALSE, stringsAsFactor=FALSE)
  return(count.vec)
})
count_matrix<-do.call(cbind, count_matrix)

##take a look on data
head(count_matrix)
tail(count_matrix)

# exclude the first 4 rows (they are aligment informaion)
count_matrix<-count_matrix[-c(1:4),]
head(count_matrix)
tail(count_matrix)
dim(count_matrix)
# [1] 8127   12


#the first column of each count.file has the gene name
#set gene names as row names
row.names(count_matrix)<-count_matrix[,1]
# keep counts only
count_matrix<-count_matrix[,c(2,4,6,8,10,12)]
#set names of count files as column names
samplenames<-gsub(".counts","",count.files)
colnames(count_matrix)<-samplenames

#explore the count matrix
head(count_matrix)
##how meny trasncript do we have?
dim(count_matrix)
# [1] 8127   6

# build 2 data frames. The first one will contain sample information and the second gene annotation
sample_data<-data.frame(sample=names(count_matrix), condition=c(rep("A", times=3), rep("B", times=3)))


# We must ensure that all genes have 10 counts in at least one condition
summary(rowSums(count_matrix[,1:3]))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0     535     874    1514    1442  447854 
summary(rowSums(count_matrix[,4:6]))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.0    348.0    554.0    982.6    882.0 262334.0


# build a index of genes with good counts
#in this case good counts are defined as those satisfying count >= 10
index_good_counts<-sapply(1:nrow(count_matrix), function(x){
  idx<-all(count_matrix[x,1:3] > rep(9, 3)) | all(count_matrix[x,4:6] > rep(9,3))
  return(idx)
})
table(index_good_counts)
#index_good_counts
#FALSE  TRUE 
#119  8008

####################
## comment on this number
# keep those genes that have good counts
count_matrix<-count_matrix[index_good_counts,]
dim(count_matrix)
#[1] 8008   6

######################
## Data exploration ##
######################
# ATENTION!!
head(sample_data$condition)
#[1] "A"  "A"  "A"  "B" "B" "B"
# Control condition must appear as a level before treatment conditions in order to compare Treat Vs. Control
sample_data$condition<-factor(sample_data$condition, levels=c("A", "B"))

# Build a DESeqDataset object
# source("http://bioconductor.org/biocLite.R")
# biocLite("DESeq2")
library(DESeq2)
library(EnhancedVolcano)

# ?DESeqDataSetFromMatrix
# We use count_matrix as countData, sample_data as colData and we choose condition column of sample_data as principal factor to the Binomial Negative Model
y_DESeq<-DESeqDataSetFromMatrix(countData=count_matrix, colData=sample_data, design= ~condition)
#set size factor equal 1 (no library size effect) to avoind rlogTransformation does for you
sizeFactors(y_DESeq)<-rep(1,6)

# To the first exploration, we make a log2 conversion over count data
#this data does not have any factor correction
matrix_rlog_DESeq<-rlogTransformation(y_DESeq, blind=TRUE)

head(count_matrix_rlog_DESeq<-assays(matrix_rlog_DESeq)[[1]])
#                La_A_1  La_A_2  La_A_3 La_B_1 La_B_2 La_B_3
#LAMA_000005000 7.238852 7.234881 7.562341      7.392334      7.151771      7.281479
#LAMA_000005100 8.003923 8.184518 8.420006      8.230631      7.981241      8.091778
#LAMA_000005200 8.884087 8.916221 9.307741      9.084360      8.908291      8.866203
#LAMA_000005300 8.388280 8.354943 8.622516      8.449911      8.261848      8.272205
#LAMA_000005400 7.742977 7.800290 8.031080      8.096424      7.805834      7.872510
#LAMA_000005500 7.541295 7.544598 7.996367      7.756523      7.511037      7.521231
#names(count_matrix_rlog_DESeq)<-names(count_matrix)

#takes a look on data distributions

###################################
##PCA over log2 transformed data ##
###################################
# The first exploration must be done over our count_matrix in order to determine if our samples are clustered togheter using count values
# Principal Component Analysis: explore samples separability  
pr_comp_y<-prcomp(t((count_matrix_rlog_DESeq)))
resumen <- summary(pr_comp_y)
labX <- signif((resumen$importance[2,1])*100, 3)
labY <- signif((resumen$importance[2,2])*100, 3)
#svg(file="PCA_DESeq_counts.svg", width=15, height=15)
par(mfrow=c(1,2))
plot(x=pr_comp_y$x[,1], y=pr_comp_y$x[,2], col=c(rep("red",3), rep("blue",3), rep("pink", 3), rep("black",3), rep("green",3)), xlab=paste("PC1 (",labX,"%)", sep=""), ylab=paste("PC2 (",labY,"%)", sep=""))
abline(h=0)
abline(v=0)
#title(main="PCA DESeq counts")
biplot(pr_comp_y,cex=c(0.5, 0.3),xlabs = count.files, ylabs = rep("", nrow(count_matrix_rlog_DESeq)), var.axes = FALSE)
#biplot(pr_comp_y,cex=c(1,0.1),xlabs=count.files,ylabs=rep("",nrow(count_matrix_rlog_DESeq)),var.axes=FALSE)
#title(main="Biplot DESeq counts")
#dev.off()
# samples are separated according the condition

#######################################
##Boxplot over log2 transformed data ##
#######################################
# The second exploration is over all counts for each sample in order to determine if samples are comparable
# Boxplot of the raw and normalized read counts
# png(file="comparisson_of_boxplot_distributions.png", width=2000, height=1200)
par(mfrow=c(1,2))
boxplot(count_matrix_rlog_DESeq, col=colors()[c(137, 134, 59,30,132,125)], main="DESeq", names=count.files)
FDP_DESeq<-sapply(1:ncol(count_matrix_rlog_DESeq), function(x){
  fdp<-density(count_matrix_rlog_DESeq[,x])
  return(fdp)})
FDP_DESeq[,1]$call<-"log2counts Densities"
plot(FDP_DESeq[,1], col=colors()[137], xlab="log2DESeqCounts", ylab="Density", xlim=c(0,16), ylim=c(0,0.4),type="l", lwd=4)
lines(FDP_DESeq[,2],col=colors()[134], lwd=4)
lines(FDP_DESeq[,3],col=colors()[59], lwd=4)
lines(FDP_DESeq[,4],col=colors()[30], lwd=4)
lines(FDP_DESeq[,5],col=colors()[132], lwd=4)
lines(FDP_DESeq[,6],col=colors()[125], lwd=4)
legend(legend=c("A1","A2","A3","B1","B2","B3"),x="topright", fill = colors()[c(137, 134, 59,30,132,125)])
# dev.off()

################################
## Library size normalization ##
################################
# There are some diferences between density functions. We should remember that existing differences in library size
# In order to estimate those diferences, we use estimateSizeFactors function
y_DESeq<-estimateSizeFactors(y_DESeq)
sizeFactors(y_DESeq)
#La_A_1   La_A_2   La_A_3    La_B_1    La_B_2    La_B_3 
#1.1064536 1.2008545 2.0391584 0.2979741 0.7520338 1.6614240 

norm_rlog_DESeqData<-rlogTransformation(y_DESeq, blind=TRUE)
count_matrix_norm_rlog_DESeqData<-assays(norm_rlog_DESeqData)[[1]]
colnames(count_matrix_norm_rlog_DESeqData)<-sample_data$sample
head(count_matrix_norm_rlog_DESeqData)


pr_comp_y<-prcomp(t((count_matrix_norm_rlog_DESeqData)))
resumen <- summary(pr_comp_y)
labX <- signif((resumen$importance[2,1])*100, 3)
labY <- signif((resumen$importance[2,2])*100, 3)
#svg(file="PCA_DESeq_counts.svg", width=15, height=15)
par(mfrow=c(1,2))
plot(x=pr_comp_y$x[,1], y=pr_comp_y$x[,2], col=c(rep("red",3), rep("blue",3), rep("pink", 3), rep("black",3), rep("green",3)), xlab=paste("PC1 (",labX,"%)", sep=""), ylab=paste("PC2 (",labY,"%)", sep=""))
abline(h=0)
abline(v=0)
#title(main="PCA DESeq counts")
biplot(pr_comp_y,cex=c(0.5, 0.3),xlabs = count.files, ylabs = rep("", nrow(count_matrix_rlog_DESeq)), var.axes = FALSE)
#biplot(pr_comp_y,cex=c(1,0.1),xlabs=count.files,ylabs=rep("",nrow(count_matrix_rlog_DESeq)),var.axes=FALSE)
#title(main="Biplot DESeq counts")
#dev.off()


#svg(file="comparisson_of_boxplot_distributions.svg", width=15, height=10)
par(mfrow=c(2,2))
boxplot(count_matrix_rlog_DESeq, col=colors()[c(137,59,101,30,125,128)], main="Rlog Counts", names=count.files)

FDP<-sapply(1:ncol(count_matrix_rlog_DESeq), function(x){
  fdp<-density(count_matrix_rlog_DESeq[,x])
  return(fdp)})

FDP_DESeq<-sapply(1:ncol(count_matrix_norm_rlog_DESeqData), function(x){
  fdp<-density(count_matrix_norm_rlog_DESeqData[,x])
  return(fdp)})

FDP[,1]$call<-"log2counts Densities"
plot(FDP[,1], col=colors()[137], xlab="log2DESeqcounts", ylab="Density", xlim=c(0,16), ylim=c(0,0.40),type="l", lwd=4)
lines(FDP[,2],col=colors()[59], lwd=4)
lines(FDP[,3],col=colors()[101], lwd=4)
lines(FDP[,4],col=colors()[30], lwd=4)
lines(FDP[,5],col=colors()[125], lwd=4)
lines(FDP[,6],col=colors()[128], lwd=4)
#legend(legend=c("A1","A2","A3","B1","B2","B3"),x="topright", fill = colors()[c(137,59,101,30,125,128)])

boxplot(count_matrix_norm_rlog_DESeqData, col=colors()[c(137,59,101,30,125,128)], main="Normalized Rlog Counts", names=count.files)

FDP_DESeq[,1]$call<-"log2counts Densities"
plot(FDP_DESeq[,1], col=colors()[137], xlab="log2DESeqcounts", ylab="Density", xlim=c(0,16), ylim=c(0,0.40),type="l", lwd=4)
lines(FDP_DESeq[,2],col=colors()[134], lwd=4)
lines(FDP_DESeq[,3],col=colors()[59], lwd=4)
lines(FDP_DESeq[,4],col=colors()[30], lwd=4)
lines(FDP_DESeq[,5],col=colors()[132], lwd=4)
lines(FDP_DESeq[,6],col=colors()[125], lwd=4)
#legend(legend=c("A1","A2","A3","B1","B2","B3"),x="topright", fill = colors()[c(137, 59,101,30,125,128)])
#dev.off()

#########################################
## Binomial Negative Model and DE Test ##
#########################################
# Before DE testing is necessary estimate dispersion of the model
y_DESeq<-estimateDispersions(y_DESeq, fitType="local")
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
par(mfrow=c(1,1))
#svg("dispersion_deseq2.svg",width=15, height=10)
plotDispEsts(y_DESeq)
#dev.off()

#Test de DE
et_DESeq<-nbinomWaldTest(y_DESeq)
results.DESeq <- results(et_DESeq)
head(results.DESeq)
#id.deseg2 <- results.DESeq$padj < 0.0001 & abs(results.DESeq$log2FoldChange) > 5
#sum(id.deseg2,na.rm=T)
head(results(et_DESeq)[order(results(et_DESeq)$padj),])
head(results(et_DESeq)[order(results(et_DESeq)$log2FoldChange),])
#log2 fold change (MLE): condition mut vs ctrl 
#Wald test p-value: condition mut vs ctrl 
#DataFrame with 6 rows and 6 columns
#baseMean log2FoldChange     lfcSE      stat       pvalue         padj
#<numeric>      <numeric> <numeric> <numeric>    <numeric>    <numeric>
#  LmxM.17.1250    3176.17        3.74551 0.1280204   29.2571 3.64421e-188 3.00757e-184
#LmxM.11.0660a   1310.27        1.65400 0.0597291   27.6917 8.79681e-169 3.63000e-165
#LmxM.34.1140    1974.53        1.30797 0.0567783   23.0365 2.01060e-117 5.53115e-114
#LmxM.30.2330    1013.39        2.74551 0.1240535   22.1317 1.56615e-108 3.23135e-105
#LmxM.04.0625    2829.12        1.48464 0.0674266   22.0186 1.91259e-107 3.15692e-104
#LmxM.26.1710    2216.84        1.15386 0.0562437   20.5153  1.57179e-93  2.16199e-90


# Outlier's analysis. For outlier detection, we should determine which of all genes has its pvalue as NA
table(is.na(results(et_DESeq)$pvalue))
table(is.na(results(et_DESeq)$padj))
# FALSE 
# 8253
#DESeq2 package allow us replace them using the trimmed mean over all samples and adjusted by the size factor
#If you found less than 223 maybe your DESeq2 version is more new and include the
#outlier detection and correction inside the nbinomWaldTest function call. 

# Following DESeq2 suggestion replace the outliers and re-fit 
#et_DESeq_clean<-replaceOutliersWithTrimmedMean(et_DESeq)
# # We must reestimate size factors, dispersion and make DE Test
#et_DESeq_clean<-estimateSizeFactors(et_DESeq_clean)
#et_DESeq_clean<-estimateDispersions(et_DESeq_clean, fitType="local")
#et_DESeq_clean<-nbinomWaldTest(et_DESeq_clean)
#check outlier genes
#table(is.na(results(et_DESeq_clean)$padj))
# FALSE  TRUE 
# 13541     8
# Note the outlier count decreasing
#a simple function to decide between up/down regulated genes
decide<-function(dds, p.adj=TRUE, value=0.05, FC_teste=1,cond1,cond2){
  if(class(dds) != "DESeqDataSet") stop("dds must be a DESeqDataSet")
  var<-"padj"
  if(!p.adj) var<-"pvalue"
  dds<-as.data.frame(results(dds,contrast = c("condition",cond1,cond2)))
  de_dds<-rep(0, nrow(dds))
  de_dds<-sapply(1:nrow(dds), function(x){
    
    if(!is.na(dds[x, var]) & dds[x, var] <= value & dds[x,"log2FoldChange"] >= FC_teste ) {
      aux<-1} else{
        if(!is.na(dds[x, var]) & dds[x, var] <= value & dds[x,"log2FoldChange"] <= -(FC_teste)) {
          aux<-(-1)
        }else {
          if(is.na(dds[x, var])) {
            aux<-NA}else aux<-0
        }}
    return(aux)
  })
  return(de_dds)
}
# We consider DE genes those that have a adjusted pavalue lower than 0.01
# table(de_DESeq <-decide(et_DESeq_clean, p.adj=TRUE, value=0.01 ))
table(de_DESeq <-decide(et_DESeq, p.adj=TRUE, value=0.05,FC_teste = 1,"RES","WT" ))
#   -1    0    1 
# 8 7709  135 
DEResults<-results(et_DESeq,contrast = c("condition","RES","WT"))
annotation<-read.table("GenesByTaxon_Summary.txt", header = TRUE, sep = "\t")
commonID<-match(row.names(DEResults),annotation$Gene_ID)
DEResults$Annotation<-annotation$Product_Description[commonID]
DEResults$Length<-annotation$Genomic_Length[commonID]


write.table(DEResults,"WT_vs_RES/DEanalysis_WTvsRES.tab", sep = "\t")
############################################
#existing_row_names <- rownames(DEResults)

#match_indices <- match(existing_row_names, geneSymbol$ID)
#new_row_names <- ifelse(is.na(match_indices), existing_row_names, geneSymbol$Gene_Symbol[match_indices])
#row.names(DEResults) <- new_row_names

#svg("volcanoPlot.svg",width=15, height=10)
EnhancedVolcano(DEResults,
                lab = rownames(DEResults),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-2,2),
                #ylim = c(0,50),
                #selectLab=c("Cqm1","Pant","Anky","CHKov1","Vnn1","CPIJ012697","CPIJ014435","NTF2","NPC2","PDhE1","39S_SRP","CPIJ002103","GST-theta","Hyp37.3-2","NA_CL_AA","CarbPepB","NOMPC","Ser3DH","CPIJ013861","39SRP","CPIJ000093","CPIJ000094","Deoxy_I","CPIJ000500"),
                pCutoff = 0.05,
                pCutoffCol = 'padj',
                FCcutoff = 1,
                pointSize = 3.0,
                labSize = 3.0)
#dev.off()
#drawConnectors = TRUE,
#max.overlaps = Inf,
#col=c('blue', 'blue', 'blue', 'red3'))#
#


############ Explore DE results ############ 
library(pheatmap)
## 1)  check if samples are clustered together using counts of DE genes
#library(gplots)
select<-which(de_DESeq ==1 | de_DESeq== -1)
pheatmap(scale(t(as.matrix(log2(counts(et_DESeq[select,], normalized=TRUE)+1))),center = TRUE,scale = TRUE))


## 2) check expression bias. 
# png(file="MA_DE_genes_DESeq.png", width=1200, height=900)
plot(rowMeans(log2(counts(et_DESeq, normalized=TRUE))),DEResults$log2FoldChange, xlab="Mean log2normalizedCounts", ylab="log2FoldChange", col=8,cex=0.5)
points(rowMeans(log2(counts(et_DESeq, normalized=TRUE)))[de_DESeq==1], DEResults$log2FoldChange[de_DESeq==1], col=3,cex=0.8)
points(rowMeans(log2(counts(et_DESeq, normalized=TRUE)))[de_DESeq==-1], DEResults$log2FoldChange[de_DESeq==-1], col=2,cex=0.8)
title(main="MA plot Log DESeq_counts vs log2FC")
legend(legend=c("UP Genes","DOWN Genes"),x="topright", fill = c(3,2))
abline(h=0, col=1, lwd=2)
# dev.off()




## 3) check length bias
# png(file="log2FCvsgenelength_distribution.png",width=1024, height=860)
plot(log2(DEResults$Length),DEResults$log2FoldChange,xlab="log2Gene_Length", ylab="log2FC", col=8, cex=0.8)
points(log2(DEResults$Length)[de_DESeq == 1],DEResults$log2FoldChange[de_DESeq== 1], col=3, cex=0.8)
points(log2(DEResults$Length)[de_DESeq == -1],DEResults$log2FoldChange[de_DESeq == -1], col=2, cex=0.8)
abline(h=0,col=1)
title("log2FC vs log2Gene_Length")
legend(legend=c("UP Genes","DOWN Genes"),x="topright", fill = c(3,2))
# dev.off()