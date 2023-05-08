library (openxlsx)

m3 <- read.xlsx ("GLONG_vs_OLA_2023_selected_mice_3_months.xlsx")
colnames (m3) <- paste (colnames (m3), "-3months-2022", sep="")

m6 <- read.xlsx ("GLONG_vs_OLA_2023_selected_mice_6_months.xlsx")
colnames (m6) <- paste (colnames (m6), "-6months-2022", sep="")

m15 <- read.xlsx ("GLONG_vs_OLA_2023_selected_mice_15_months.xlsx")
colnames (m15) <- paste (colnames (m15), "-15months-2022", sep="")


anno <- m3[ ,c("gene_name-3months-2022", "gene_type-3months-2022", "mgi_id-3months-2022", "description-3months-2022")]
colnames (anno) <- gsub ("-3months-2022", "", colnames (anno))
head (anno)



mco <- merge (m3, m6, by.x="Geneid-3months-2022", by.y="Geneid-6months-2022")
mco <- merge (mco, m15, by.x="Geneid-3months-2022", by.y="Geneid-15months-2022")

mco <- mco[ ,grep ("type|external|description|lfcSE|stat|pvalue|mgi|GLONG|OLA" , colnames (mco), invert=TRUE)]
mco <- mco [order (mco$`padj-3months-2022`), ]
colnames (mco)[colnames (mco) == "Geneid-3months-2022"] <- "Geneid"
mco$stats <- ifelse (mco$`padj-3months-2022` < 0.05 & mco$`padj-15months-2022` < 0.05 & mco$`padj-6months-2022` < 0.05, "Yes","No")
mco$trend <- ifelse (mco$`log2FoldChange-3months-2022` < 0 & mco$`log2FoldChange-15months-2022` < 0 & mco$`log2FoldChange-6months-2022` < 0, "Down","No")
mco$trend [mco$`log2FoldChange-3months-2022` > 0 & mco$`log2FoldChange-15months-2022` > 0 & mco$`log2FoldChange-6months-2022` > 0] <- "Up"
mco <- merge (mco, anno, by.x="gene_name-3months-2022", by.y="gene_name", all.x=TRUE)
mco <- mco[order (mco$`padj-3months-2022`), ]
head (mco)

write.xlsx (mco, "Comparison of GLONG vs OLA at 3, 6 and 15 months.xlsx", rowNames=F)


