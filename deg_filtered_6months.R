## 220120_A00558_0160_AH5V77DSX3

library (openxlsx)
library (DESeq2)
library (ggplot2)
library (ggrepel)
library (pheatmap)


anno <- read.delim ("gencode.vM32.annotation.txt")
anno <- anno[ ,grep ("transcript_id", colnames (anno), invert=TRUE)]
anno <- unique (anno)

a <- read.delim ("subread.counts.txt", skip=1)
a <- a[ ,grep ("ene|bam", colnames (a))]
a <- a[grep ("miRNA|Mt_tRNA|Mt_rRNA|rRNA|snRNA|snoRNA|scRNA|sRNA|misc_RNA|scaRNA|ribozyme|IG_|TR_", a$gene_type, invert=TRUE), ]
colnames (a) <- gsub ("_S[0-9]+.*", "", colnames (a)) 
colnames (a) <- gsub (".*TRM_", "", colnames (a))

a <- merge (a, anno, by.x="Geneid", by.y="gene_id", all.x=TRUE) 
a <- a[ ,grep ("gene_type.y|gene_name.y", colnames (a), invert=TRUE)]
colnames (a) [colnames (a) == "gene_name.x"] <- "gene_name"
colnames (a) [colnames (a) == "gene_type.x"] <- "gene_type"

#write.xlsx (a, "star_gene_raw_counts.xlsx", rowNames=F)

annot <- a
annot <- annot[ ,c("Geneid", "gene_name", "gene_type", "mgi_id", "external_gene_name", "description")]


torm <- c("gene_name", "gene_type", "mgi_id", "external_gene_name", "description")
a <- a[ ,!colnames (a) %in% torm]
row.names (a) <- a[ ,1]
colnames (a) <- gsub ("star.IIT_RNM_", "", colnames (a))
a <- a[ ,-1]


pheno <- read.delim ("Moratalla phenodata 220303.txt")
pheno <- pheno[pheno$area == "SN", ]
pheno <- pheno[pheno$experiment ==2, ]
idx <- match (colnames (a), pheno$sample)
pheno <- pheno[idx, ]


stopifnot (colnames (a) == pheno$sample)

# 6 months old
age <- 6

pheno.s <- pheno[pheno$age == age, ]

pheno.s <- pheno[pheno$sample %in% c("GLONG13SN", "GLONG17SN", "GLONG9SN", "OLA1SN", "OLA2SN"), ]

a.s <- a[ ,colnames (a) %in% pheno.s$sample]
stopifnot (colnames (a.s) == pheno.s$sample)

dds <- DESeqDataSetFromMatrix(countData = round (a.s), colData = pheno.s, design = ~ genotype)

# Three samples only !!
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep,]
dds

dds <- DESeq(dds)
res <- results(dds, contrast=c("genotype", "GLONG", "OLA"))

res <- merge (data.frame (res), counts (dds), by="row.names")
#res <- merge (data.frame (res), round (counts (dds, normalized=TRUE)), by="row.names")
res <- merge (res, annot, by.x="Row.names", by.y="Geneid")
colnames (res)[1] <- "Geneid"
res <- res[order (res$padj), ]

# Sanity check
res[res$gene_name == "Snca", ] 
res[res$gene_name == "Th", ] 

table (res$padj < 0.05)

write.xlsx (res, paste (paste ("GLONG_vs_OLA_2023_selected_mice_" , age, sep=""), "_months.xlsx", sep=""), rowNames=F)


## heatmap plot
# Th, DDC, NR4A2, SLC6A3, SNCA, Hmgn2, Pafah1b2, Srp54A, NUDC-PS1
select <- c("ENSMUSG00000000214.12", "ENSMUSG00000020182.17", "ENSMUSG00000026826.14", "ENSMUSG00000021609.7", "ENSMUSG00000025889.14", "ENSMUSG00000003038.16", "ENSMUSG00000003131.8", "ENSMUSG00000073079.7", "ENSMUSG00000110331.2")

df <- as.data.frame(colData(dds)[,c("genotype","sample")])

pdf (paste (paste ("Heatmap selected plot", age, sep=" "), "months.pdf", sep=" "))
pheatmap( log2 (counts(dds,normalized=TRUE)+1) [row.names (counts(dds)) %in% select, ], cluster_rows=FALSE, show_rownames=TRUE, cluster_cols=FALSE, annotation_col=df)
dev.off ()


## PCA plot
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("genotype", "sample"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=genotype, label=sample)) +
  		geom_point(size=3) +
  		xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  		ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  		geom_text_repel()  + 
		  coord_fixed () 

ggsave (paste (paste ("PCA selected plot", age, sep=" "), "months.pdf", sep=" "))














