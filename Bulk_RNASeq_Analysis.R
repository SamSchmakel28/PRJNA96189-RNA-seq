# PRJNA96189
# https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA961819&o=acc_s%3Aa

library(DESeq2)
library(ggplot2)
library(biomaRt)

cts <- read.delim2("C:/Users/artsc/OneDrive/Documents/RNASeq_Project/PRJNA96189/PRJNA96189cts.txt", header = TRUE, skip = 1, sep = "\t")
colnames(cts)[7:16] <- c("SRR24302768","SRR24302769","SRR24302770","SRR24302771","SRR24302772",
                         "SRR24302773","SRR24302774","SRR24302775","SRR24302776","SRR24302777")

# Construct coldata table
condition <- c("Gata3_low","Gata3_low","Gata3_low","Gata3_low","Gata3_low","Gata3_high","Gata3_high",
               "Gata3_high","Gata3_high","Gata3_high")
coldata <- as.data.frame(condition)
rownames(coldata) <- colnames(cts)[7:16]

cts_only <- cts[,7:16]
rownames(cts_only) <- cts$Geneid
dds <- DESeqDataSetFromMatrix(countData = cts_only,
                              colData = coldata,
                              design = ~ condition)

smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

dds$condition <- factor(dds$condition, levels = c("Gata3_low","Gata3_high"))

dds <- DESeq(dds)
res <- results(dds)
res

# Creating normalized counts and create a PCA plot
normalizedcts <- counts(dds, normalized = TRUE)
rld <- rlogTransformation(dds)
pca <- plotPCA(rld)

sum(res$padj < 0.05, na.rm=TRUE)

#resLFC computation
#resLFC <- lfcShrink(dds, coef="condition_treated_vs_untreated", type="apeglm")
#resLFC

res05 <- results(dds, alpha=0.05)
summary(res05)

res_na_padj <- res[is.na(res$padj),]
res <- res[!is.na(res$padj),]
res <- res[res$padj < 0.05,]
res_sig_up <- res[res$log2FoldChange > 0,]
res_sig_down <- res[res$log2FoldChange < 0,]

res_sig_up_table <- data.frame(res_sig_up)
res_sig_down_table <- data.frame(res_sig_down)

# Sort by log fold change and save just the gene list
res_sig_up

#res_test <- res[is.na(res$padj),]

# Convert gene names

ensembl <- useMart('ensembl', dataset = 'mmusculus_gene_ensembl')
