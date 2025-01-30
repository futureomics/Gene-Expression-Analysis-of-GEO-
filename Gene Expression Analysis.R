# Gene Expression Analysis

if(!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GEOquery", "limma", "clusterProfiler", "org.Hs.eg.db", "enrichplot", "pheatmap", "RColorBrewer", "rmarkdown"))

library(GEOquery) 
gds <- getGEO("GDS2778") # Download GDS2778 data from GEO
eset <- GDS2eSet(gds)
pData(eset)[ ,1:2] 

dim(exprs(eset))
exprs(eset)[1:4, ] 
#DEG analysis with Limma
library(limma) # Loads limma library.
targets <- pData(eset)
design <- model.matrix(~ -1+factor(as.character(targets$agent)))
colnames(design) <- c("Treatment", "Control") 
fit <- lmFit(eset, design) 
contrast.matrix <- makeContrasts(Control-Treatment, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix) 
fit2 <- eBayes(fit2)
deg_df <- topTable(fit2, coef=1, adjust="fdr", sort.by="B", number=Inf)
deg_df <- deg_df[!is.na(deg_df$adj.P.Val),] # Removes rows with NA values
head(deg_df)[,22:27]

#DEG filtering
nr <- nrow(deg_df)
nr

pf <- nrow(deg_df[deg_df$adj.P.Val <= 0.01, ])
pf
#The number of DEGs passing this filter is 26.

fcf <- nrow(deg_df[deg_df$logFC >= 1 | deg_df$logFC <= -1, ])
fcf
#The number of DEGs passing this filter is 176.

pf_log <- deg_df$adj.P.Val <= 0.01
fcf_log <- deg_df$logFC >= 1 | deg_df$logFC <= -1
comb_filter <- pf_log & fcf_log
combf <- nrow(deg_df[comb_filter, ])
combf
#The number of DEGs passing this filter is 16.

#Enrichment analysis of GO terms
#The following performs over-representation analysis (ORA) of GO terms using functions from the clusterProfiler package (Yu et al. 2012).

#Prepare input
#For the GO term enrichment analysis the DEGs are filtered using a adjusted p-value of ≤0.05
#Subsequently, the gene identifiers (Entrez IDs) are stored in a character vector named ids

cutoff <- 0.05 # Cutoff to use for filtering on adjusted p-value (FDR)
ids <- deg_df[deg_df$adj.P.Val <= cutoff, "Gene.ID"]
ids <- ids[!grepl("/", ids)] # Removes ambiguous mappings
ids <- ids[nchar(ids)!=0] # Removes empty slots
ids # Prints gene IDs

#Enrichment of MF terms
library(clusterProfiler); library(org.Hs.eg.db); library(enrichplot)
ego_mf <- enrichGO(gene=ids, OrgDb=org.Hs.eg.db, ont="MF", pAdjustMethod="BH", pvalueCutoff=0.1, readable=TRUE)
dim(ego_mf) # Returns number of rows and columns in result table

head(ego_mf) # Returns first six rows to inspect results

#Visualization of MF result
#The barplot plots the top scoring GO terms (here 10) in form of a bar plot. 
#To plot the tree structure of the corresponding DAG, the goplot function can be used.
barplot(ego_mf, showCategory=10)

goplot(ego_mf)

#Enrichment of BP terms
#Same as above but for Biological Process (BP) Gene Ontology
ego_bp <- enrichGO(gene=ids, OrgDb=org.Hs.eg.db, ont="BP", pAdjustMethod="BH", pvalueCutoff=0.1, readable=TRUE)
dim(ego_bp)
head(ego_bp)

#Visualization of BP result
#Same bar plot as above but for Biological Process (BP) Gene Ontology
barplot(ego_bp, showCategory=10)

#Clustering
#The following uses DEGs passing an adjusted p-value cutoff of ≤0.05
#to subset the gene expression matrix imported from GEO.
library(pheatmap); library("RColorBrewer")
cutoff <- 0.05 # Cutoff to use for filtering on adjusted p-value (FDR)
affy_ids <- row.names(deg_df[deg_df$adj.P.Val <= cutoff, ])
deg_ma <- exprs(eset)[affy_ids, ] 
pheatmap(deg_ma, scale="row", color=brewer.pal(9, "Blues"))
