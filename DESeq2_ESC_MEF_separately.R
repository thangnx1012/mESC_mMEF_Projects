## RNA-seq analysis with DESeq2
setwd("/media/hkh/8TB/XUANTHANG/INO80_ESC_MEF/RNAseq/Results")

# RNA-seq data from GSE95633 INO80_MEF/GSE49252 INO80 ChIPseq_ESC/GSE158343 INO80_RNAseq_ESC
# All patients with 4 Control ESC (should be remove 2) , 2 KD_INO80 ESC, 2 control MEF and 2 KD_INO80_MEF


#1 Import counting read & pre-process ----------------------------------------------------

# Import data from featureCounts
## Previously ran at command line something like this:
## featureCounts -a genes.gtf -o INO80_KD_ESC_MEF_counts.tsv -t -T 6 -t exon -g gene_id SRR*.bam
## read count file using column names as first row and row names as first column
countdata <- read.table("featureCounts/INO80_KD_ESC_MEF_counts.tsv", header=TRUE, row.names=1) #sep="\t"
head(countdata)
# Remove first five columns (chr, start, end, strand, length)
countdata <- countdata[ ,6:ncol(countdata)]

# Remove .bam or .sam from filenames
colnames(countdata) <- gsub("\\Aligned.sortedByCoord.out.bam", "", colnames(countdata))
colnames(countdata) <- gsub("\\Results.STAR.", "", colnames(countdata))
head(countdata)
# remove 2 control from ESC
countdata <- countdata[,-c(1:2)]
countdata <- countdata[,c(1:5,7,6,8)] # rearrange column data

ESC_countdata <- countdata[,c(1:4)]
MEF_countdata <- countdata[,c(5:8)]
# Convert to matrix
countdata <- as.matrix(countdata)
head(countdata)

# Assign condition (first two are ESC controls, second one contain the expansion of each ESC_INO80KD, MEF_Cont and MEF_INO80_KD)
ESC_celltype <- factor(rep("ESC",4))
MEF_celltype <- factor(rep("MEF",4))
ESC_celltype
MEF_celltype
condition <- factor(c(rep("Cont", 2), rep("INO80_KD", 2)))

ESC_condition
MEF_condition
#2. Analysis with DESeq2 ----------------------------------------------------
library(DESeq2)

# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
ESC_coldata <- data.frame(row.names=colnames(ESC_countdata), ESC_celltype, condition)
MEF_coldata <- data.frame(row.names=colnames(MEF_countdata), ESC_celltype, condition)
ESC_coldata
MEF_coldata
ESC_dds <- DESeqDataSetFromMatrix(countData=ESC_countdata, colData=ESC_coldata, design=~condition)
MEF_dds <- DESeqDataSetFromMatrix(countData=MEF_countdata, colData=MEF_coldata, design=~condition)
ESC_dds
MEF_dds

# Pre-filtering 
keep <- rowSums(counts(ESC_dds) >= 10) >= 3 #low count gene at least 3 samples with a count of 10 or higher
ESC_dds <- ESC_dds[keep,]
keep <- rowSums(counts(MEF_dds) >= 10) >= 3 #low count gene at least 3 samples with a count of 10 or higher
MEF_dds <- MEF_dds[keep,]
rm(keep)
# Differential expression analysis

# Run the DESeq pipeline
ESC_dds <- DESeq(ESC_dds)
MEF_dds <- DESeq(MEF_dds)

# Plot dispersions
plotDispEsts(ESC_dds, main="Dispersion plot")
plotDispEsts(MEF_dds, main="Dispersion plot")
#save as png file
tiff("Rplot/qc-dispersions.tiff", 1000, 1000, pointsize=20)
plotDispEsts(ESC_dds, main="Dispersion plot")
dev.off()

# Principal components analysis
# Regularized log transformation for clustering or heatmaps, etc
rld <- rlogTransformation(dds)
head(assay(rld))
hist(assay(rld))

## Could do with built-in DESeq2 function:
## DESeq2::plotPCA(rld, intgroup=c("celltype", "condition"))


# Get differential expression results
resultsNames(ESC_dds)
res_ESC <- results(ESC_dds, contrast=c("condition","INO80_KD","Cont"))
FC = 1.5
sum(res_ESC$padj < 0.05 & res_ESC$log2FoldChange >= log2(FC), na.rm=TRUE)
sum(res_ESC$padj < 0.05 & res_ESC$log2FoldChange <= -log2(FC), na.rm=TRUE)

## Order by adjusted p-value
resOrdered_ESC <- res_ESC[order(res_ESC$padj),]
## Extract DEG_all
sig <- resOrdered_ESC[!is.na(resOrdered_ESC$padj) & resOrdered_ESC$padj < 0.05 &
                    abs(resOrdered_ESC$log2FoldChange)>=log2(FC),]
ESC_DEG_All <- rownames(sig)
## extract DEG_up_regulated gene list
up <- resOrdered_ESC[!is.na(resOrdered_ESC$padj) & resOrdered_ESC$padj < 0.05 &
                       resOrdered_ESC$log2FoldChange >= log2(FC),]
ESC_DEG_up_regulated_gene <- rownames(up)

## extract DEG_down_regulated gene list
down <- resOrdered_ESC[!is.na(resOrdered_ESC$padj) & resOrdered_ESC$padj < 0.05 &
                         resOrdered_ESC$log2FoldChange <= -log2(FC),]
ESC_DEG_down_regulated_gene <- rownames(down)




resultsNames(MEF_dds)
res_MEF <- results(MEF_dds, contrast=c("condition","INO80_KD","Cont"))
FC = 1.5
sum(res_MEF$padj < 0.05 & res_MEF$log2FoldChange >= log2(FC), na.rm=TRUE)
sum(res_MEF$padj < 0.05 & res_MEF$log2FoldChange <= -log2(FC), na.rm=TRUE)

## Order by adjusted p-value
resOrdered_MEF <- res_MEF[order(res_MEF$padj),]
## Extract DEG_all
sig <- resOrdered_MEF[!is.na(resOrdered_MEF$padj) & resOrdered_MEF$padj < 0.05 &
                        abs(resOrdered_MEF$log2FoldChange)>=log2(FC),]
MEF_DEG_All <- rownames(sig)
## extract DEG_up_regulated gene list
up <- resOrdered_MEF[!is.na(resOrdered_MEF$padj) & resOrdered_MEF$padj < 0.05 &
                       resOrdered_MEF$log2FoldChange >= log2(FC),]
MEF_DEG_up_regulated_gene <- rownames(up)

## extract DEG_down_regulated gene list
down <- resOrdered_MEF[!is.na(resOrdered_MEF$padj) & resOrdered_MEF$padj < 0.05 &
                         resOrdered_MEF$log2FoldChange <= -log2(FC),]
MEF_DEG_down_regulated_gene <- rownames(down)

All_DEG <- intersect(ESC_DEG_All, MEF_DEG_All)

All_down <- intersect(ESC_DEG_down_regulated_gene, MEF_DEG_down_regulated_gene)

All_Up <- intersect(ESC_DEG_up_regulated_gene, MEF_DEG_up_regulated_gene)








write.csv(DEG_All, file="DEG_All_results.csv", quote = F, row.names = F)














## Examine plot of p-values
hist(res$pvalue, breaks=50, col="grey")

## Examine independent filtering
attr(res, "filterThreshold")
plot(attr(res,"filterNumRej"), type="b", xlab="quantiles of baseMean", ylab="number of rejections")























# Colors for plots below
## Ugly:
## (mycols <- 1:length(unique(condition)))
## Use RColorBrewer, better
library(RColorBrewer)
mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition))]
# Sample distance heatmap
sampleDists <- as.matrix(dist(t(assay(rld))))
library(gplots)
png("RPlot/qc-heatmap-samples.png", w=1000, h=1000, pointsize=20)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "Blue", "white"),
          ColSideColors=mycols[condition], RowSideColors=mycols[condition],
          margin=c(10, 10), main="Sample Distance Matrix")
dev.off()

df <- as.data.frame(colData(dds)[,c("condition","celltype")])


#pheatmap(normalized_counts[DEG_All,], 
#             color=greenred(75),
#             fontsize_row=1,scale="row",
#             cluster_rows=TRUE, show_rownames=TRUE,
#             cluster_cols=TRUE, annotation_col=df, )


vsd <- vst(dds, blind=FALSE)
#using vst normalization(recommended)
library(pheatmap)
pl<-pheatmap(assay(vsd)[DEG_All,], 
             color=greenred(75),
             scale="row",
             fontsize_row=1, cutree_rows = 2,
             cluster_rows=TRUE, show_rownames=TRUE, 
             cluster_cols=FALSE, annotation_col=df, )


hc <- pl$tree_row
lbl <- cutree(hc, 2)
cluster1<-which(lbl==1)
cluster2<-which(lbl==2)

write.table(cluster1, "cluster1.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(cluster2, "cluster2.txt", sep = "\t", quote = F, row.names = F, col.names = F)



## MA plot
## Merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)
## Could do with built-in DESeq2 function:
## DESeq2::plotMA(dds, ylim=c(-1,1), cex=1)
## I like mine better:
maplot <- function (res, thresh=0.05, labelsig=TRUE, textcx=1, ...) {
  with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
  with(subset(res, padj<thresh), points(baseMean, log2FoldChange, col="red", pch=20, cex=1.5))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<thresh), textxy(baseMean, log2FoldChange, labs=Gene, cex=textcx, col=2))
  }
}
png("RPlot/diffexpr-maplot.png", 1500, 1000, pointsize=20)
maplot(resdata, main="MA Plot")
dev.off()

## Volcano plot with "significant" genes labeled
volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}
png("RPlot/diffexpr-volcanoplot.png", 1200, 1000, pointsize=20)
volcanoplot(resdata, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-2.3, 2))
dev.off()

































## I like mine better:
rld_pca <- function (rld, intgroup = "condition", ntop = 500, colors=NULL, legendpos="bottomleft", main="PCA Biplot", textcx=1, ...) {
  require(genefilter)
  require(calibrate)
  require(RColorBrewer)
  rv = rowVars(assay(rld))
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld)[select, ]))
  fac = factor(apply(as.data.frame(colData(rld)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
  if (is.null(colors)) {
    if (nlevels(fac) >= 3) {
      colors = brewer.pal(nlevels(fac), "Paired")
    }   else {
      colors = c("black", "red")
    }
  }
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC1 (",as.character(pc2var),"%)")
  plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)
  with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
  legend(legendpos, legend=levels(fac), col=colors, pch=20)
  #     rldyplot(PC2 ~ PC1, groups = fac, data = as.data.frame(pca$rld),
  #            pch = 16, cerld = 2, aspect = "iso", col = colours, main = draw.key(key = list(rect = list(col = colours),
  #                                                                                         terldt = list(levels(fac)), rep = FALSE)))
}
png("qc-pca.png", 1000, 1000, pointsize=20)
rld_pca(rld, colors=mycols, intgroup="condition", xlim=c(-75, 35))
dev.off()
