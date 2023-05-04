library (openxlsx)

m3 <- read.xlsx ("GLONG_vs_OLA_2023_selected_mice_3_months.xlsx")
colnames (m3) <- paste (colnames (m3), "-3months-2022", sep="")

m15 <- read.xlsx ("GLONG_vs_OLA_2023_selected_mice_15_months.xlsx")
colnames (m15) <- paste (colnames (m15), "-15months-2022", sep="")

mco <- merge (m3, m15, by.x="Geneid-3months-2022", by.y="Geneid-15months-2022", all.x=TRUE, all.y=TRUE)
mco <- mco[ ,grep ("type|external|description|lfcSE|stat|pvalue|mgi|GLONG|OLA" , colnames (mco), invert=TRUE)]
mco <- mco [order (mco$`padj-3months-2022`), ]
colnames (mco)[colnames (mco) == "Geneid-3months-2022"] <- "Geneid"
mco$stats <- ifelse (mco$`padj-3months-2022` < 0.05 & mco$`padj-15months-2022` < 0.05, "Yes","No")
mco$trend <- ifelse (mco$`log2FoldChange-3months-2022` < 0.05 & mco$`log2FoldChange-15months-2022` < 0.05, "Yes","No")
write.xlsx (mco, "Comparison of GLONG vs OLA at 3 and 15 months.xlsx", rowNames=F)


