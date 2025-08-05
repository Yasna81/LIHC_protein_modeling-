#library(jsonlite)
#json_data <- fromJSON("data/TCGA.LIHC.sampleMap_HiSeqV2.json")
#str(json_data)

exp <- read.delim("data/TCGA.LIHC.sampleMap_HiSeqV2",
                  header = TRUE,
                  row.names = 1 ,
                  check.names = FALSE)
dim(exp)
sample_ids <- colnames(exp)
sample_type_code <- substr(sample_ids,14,15)
#we assign labels
sample_group <- ifelse(sample_type_code == "01", "Tumor",
                       ifelse(sample_type_code == "11","Normal", NA))
metadata <- data.frame(
    sample_id = sample_ids,
    condition = sample_group
)

metadata <- na.omit(metadata)
exp_matrix_filtered <- exp[,metadata$sample_id]
#now metadata and expr matrix is ready 
#BiocManager::install("DESeq2")
library(DESeq2)
metadata <- metadata[match(colnames(exp_matrix_filtered),metadata$sample_id),]
count_data <- round(as.matrix(exp_matrix_filtered))
dds <- DESeqDataSetFromMatrix(
    countData =  count_data,
    colData = metadata ,
    design = ~ condition 
)

dds <- dds[rowSums(counts(dds)) > 10 ,]
dds <- DESeq(dds)
res <- results(dds , contrast =c("condition" , "Tumor" , "Normal"))
deg_df <- as.data.frame(res)
deg_df <- deg_df[!is.na(deg_df$padj),]
degs <- rownames(deg_df[deg_df$padj < 0.05 & abs(deg_df$log2FoldChange) > 1 , ])

#filtering expression 
exp_degs <- exp_matrix_filtered[rownames(exp_matrix_filtered) %in% degs,]
#pr pr interactions 
library(STRINGdb)
string_db <- STRINGdb$new(version = "11.5", species = 9606 , score_threshold = 700)
degs_data  <- data.frame(gene = degs)
mapped_degs <- string_db$map(degs_data, "gene", removeUnmappedRows = TRUE)
#getting interactions in degs 
ppi_data <- string_db$get_interactions(mapped_degs$STRING_id)
# we did that manually 
ppi_data <- read.table("data/9606.protein.links.v11.5.txt",header = TRUE,sep = "", stringsAsFactors = FALSE)
head(ppi_data)
ppi_data <- subset(ppi_data,combined_score >= 700)
valid_ids <- mapped_degs$STRING_id
ppi_filtered <- subset (ppi_data , protein1 %in% valid_ids & protein2 %in% valid_ids)
#no net 
id_to_gene <- mapped_degs[,c("gene","STRING_id")]
ppi_merged <- merge(ppi_filtered, id_to_gene, by.x="protein1",by.y = "STRING_id")
ppi_merged <-merge(ppi_merged, id_to_gene,by.x = "protein1" , by.y= "STRING_id", suffixes = c("_gene1","_gene2"))
pos_pair <- data.frame(gene1 = ppi_merged$gene_gene1,
                       gene2 = ppi_merged$gene_gene2,
                       label = 1)

#corrolation
#for better we need to standerdized 
vsd <- vst(dds,blind = TRUE)
vsd_mat <- assay(vsd)
deg_genes <- rownames(exp_degs)
vsd_deg <- vsd_mat[deg_genes,]
vsd_deg_scaled <- t(scale(t(vsd_deg)))

#ppi to gene ids 
library(biomaRt)
ensembl <- useMart("ensembl",dataset = "hsapiens_gene_ensembl")
ppi_data$protein1_clean <- sub("9606\\.","",ppi_data$protein1)
ppi_data$protein2_clean <- sub("9606\\.","",ppi_data$protein2)
all_prot <- unique(c(ppi_data$protein1_clean, ppi_data$protein2_clean))
prot_map <- getBM(attributes = c("ensembl_peptide_id","hgnc_symbol","entrezgene_id"),
                  filters = "ensembl_peptide_id",
                  values = all_prot,
                  mart = ensembl)
#1
ppi_data <- merge(ppi_data,prot_map, by.x = "protein1_clean", by.y = "ensembl_peptide_id",all.x = TRUE)
names(ppi_data)[names(ppi_data)== "hgnc_symbol"] <- "gene1"
names(ppi_data)[names(ppi_data) == "entrezgene_id"] <- "entrez1"
#2
ppi_data <- merge(ppi_data,prot_map, by.x = "protein2_clean", by.y = "ensembl_peptide_id",all.x = TRUE)
names(ppi_data)[names(ppi_data)== "hgnc_symbol"] <- "gene2"
names(ppi_data)[names(ppi_data) == "entrezgene_id"] <- "entrez2"
#we have 517000 gene pairs, to have it more manageble we randomly choose 5000 with combined score more than 900
set.seed(42)
sample_ppi <- ppi_data[sample(nrow(ppi_data), 5000),]
filterd_ppi <- subset(sample_ppi, combined_score > 900)

#GGO semantic 
library(GOSemSim)
library(org.Hs.eg.db)
hsGO <- godata('org.Hs.eg.db',ont = "BP",computeIC = TRUE)
filterd_ppi$entrez1 <- as.character(filterd_ppi$entrez1)
filterd_ppi$entrez2 <- as.character(filterd_ppi$entrez2)
go_sim <- vector("numeric", length = nrow(filterd_ppi))

for (i in 1:nrow(filterd_ppi)) {
    id1 <- filterd_ppi$entrez1[i]
    id2 <- filterd_ppi$entrez2[i]
    
    if (!is.na(id1) && !is.na(id2)) {
        sim_result <- tryCatch({
            geneSim(id1, id2, semData = hsGO, measure = "Wang")
        }, error = function(e) return(NA))
        
        if (is.list(sim_result) && !is.null(sim_result$geneSim)) {
            go_sim[i] <- sim_result$geneSim
        } else if (is.numeric(sim_result)) {
            go_sim[i] <- sim_result
        } else {
            go_sim[i] <- NA
        }
    } else {
        go_sim[i] <- NA
    }
}


filterd_ppi$go_similarity <- go_sim
#labeling and correlation 
library(dplyr)
filterd_ppi_2 <- filterd_ppi %>% select(-label)
filterd_ppi_2 <- filterd_ppi_2 %>% filter(go_similarity >= 0.6 | go_similarity < 0.4) %>%  mutate(label= ifelse(go_similarity >= 0.6,1,0))
table(filterd_ppi_2$label)




set.seed(123)
train_idx <- sample(1 : nrow(filterd_ppi_2), size = 0.7 * nrow(filterd_ppi_2))
train_df <- filterd_ppi_2[train_idx,]
test_df <- filterd_ppi_2[-train_idx,]


get_correlation <- function(df, vsd_deg_scaled) {
    apply(df, 1, function(row) {
        g1 <- row["gene1"]
        g2 <- row["gene2"]
        if (g1 %in% rownames(vsd_deg_scaled) && g2 %in% rownames(vsd_deg_scaled)) {
            return(cor(as.numeric(vsd_deg_scaled[g1, ]), as.numeric(vsd_deg_scaled[g2, ])))
        } else {
            return(NA)
        }
    })
}

# Add correlation to train and test separately
train_df$correlation <- get_correlation(train_df,  exp)
test_df$correlation <- get_correlation(test_df,  exp)
#impute na corrolations 
train_df$correlation[is.na(train_df$correlation)] <- 0
test_df <- test_df[!is.na(test_df$correlation),]

model<- glm(
    label ~ correlation + combined_score,
    data = train_df,
    family = "binomial"
)
summary(model)



# Predict and evaluate
pred <- predict(model, test_df, type = "response")
library(pROC)
roc_obj <- roc(test_df$label, pred)
auc(roc_obj)
plot(roc_obj, main = paste("AUC=",round(auc(roc_obj),2)))
library(caret)
confusionMatrix(as.factor(ifelse(pred > 0.5,1,0)), as.factor(test_df$label))
#evaluate bio
test_df$pred_prob <- pred
top_pairs <- test_df[order(-test_df$pred_prob),][1:100,]
top_entrez <- unique(c(top_pairs$entrez1, top_pairs$entrez2))
library(clusterProfiler)
library(org.Hs.eg.db)
kegg_results <- enrichKEGG(gene = top_entrez,
                           organism = "hsa",
                           pvalueCutoff = 0.05)
dotplot(kegg_results, showCategory = 10)
go_reults  <- enrichGO(gene = top_entrez,
                       OrgDb = org.Hs.eg.db,
                       pvalueCutoff = 0.05,
                       ont = "BP",
                       readable = TRUE)
dotplot(go_reults, showCategory = 10)

#plots

top_genes <- deg_df[deg_df$padj < 0.05 & abs(deg_df$log2FoldChange) > 1 , ]
 
top_genes <- top_genes[order(top_genes$padj), ][1:10, ]
deg_df$significant <- "Not Significant"  
deg_df$significant[deg_df$padj < 0.05 & abs(deg_df$log2FoldChange) > 1] <- "Significant"


#network diagram 
edge_list <- data.frame(
    source = test_df$gene1,
    target = test_df$gene2,
    score = test_df$pred_prob
)
edge_list_filterd <- subset(edge_list, score > 0.8)
library(igraph)
g <- graph_from_data_frame(edge_list_filterd, directed = FALSE)
plot(g, vertex.size = 20, vertex.label.cex = 1, edge.width = edge_list_filterd$score * 7)
library(DESeq2)
vsd <- vst(dds,blind= FALSE)
vst_matrix <- assay(vsd)

# degs
top_degs <- res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1), ]
top_degs <- top_degs[order(abs(top_degs$log2FoldChange), decreasing = TRUE), ]
top_genes <- rownames(top_degs)[1:50]

heatmap_data <- vst_matrix[top_genes, ]
heatmap_data <- t(scale(t(heatmap_data)))
library(pheatmap)
pheatmap(heatmap_data,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 10)
head(top_genes)
#looking for top down and up-regulated
library(dplyr)  
library(DESeq2)
  
res <- results(dds, contrast = c("condition", "Tumor", "Normal")) %>%  
    as.data.frame() %>%  
    filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%  
    arrange(desc(log2FoldChange))   

# Top 10 upregulated genes 
top_upregulated <- res %>%  
    slice_max(log2FoldChange, n = 10)  

# Top 10 downregulated genes 
top_downregulated <- res %>%  
    slice_min(log2FoldChange, n = 10)  





#enhanced volcano
library(ggplot2)
library(ggrepel) 

# Prepare data
degs <- res %>% 
    tibble::rownames_to_column("gene") %>%  # Convert row names to column
    mutate(
        significance = case_when(
            padj < 0.05 & abs(log2FoldChange) > 1 ~ "Significant",
            TRUE ~ "Not significant"
        )
    )

# Define  label
top_genes <- rbind(
    degs %>% arrange(desc(log2FoldChange)) %>% head(10),  # Top upregulated
    degs %>% arrange(log2FoldChange) %>% head(10)         # Top downregulated
)


ggplot(degs, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
   
    geom_point(aes(color = ifelse(log2FoldChange > 0, "Up", "Down")),alpha = 0.6, size = 2) +
    
    
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
   
    geom_text_repel(
        data = top_genes,
        aes(label = gene),
        color = "black",
        size = 3,
        box.padding = 0.5,
        max.overlaps = Inf
    ) +
    
    
    scale_color_manual(values = c("blue", "red")) +
    labs(
        x = "log2 Fold Change",
        y = "-log10(Adjusted p-value)",
        color = "Significance",
        title = "Volcano Plot of DEGs"
    ) +
    theme_minimal() +
    theme(
        legend.position = "top",
        panel.grid.major = element_line(color = "gray90")
    )

#network pf prob
 

library(dplyr)
#
top_10_ppis <- test_df %>%
    arrange(desc(pred_prob)) %>%  # Sort by prediction score (highest first)
    head(20) 

library(visNetwork)

# Nodes: Unique genes in top 20 PPIs
nodes <- data.frame(
    id = unique(c(top_10_ppis$gene1, top_10_ppis$gene2)),
    label = unique(c(top_10_ppis$gene1, top_10_ppis$gene2))
)


edges <- top_10_ppis %>%
    mutate(
        from = gene1,
        to = gene2,
        width = pred_prob * 5,  
        title = paste("pred_prob:", round(pred_prob, 3))  # Tooltip
    )

visNetwork(nodes, edges) %>%
    visNodes(shape = "dot", size = 20, font = list(size = 20)) %>%
    visEdges(smooth = FALSE, arrows = "none") %>%
    visLayout(randomSeed = 42)  # Consistent layout


#github saving 

write.csv(test_df, file = "~/project/output_data.csv", row.names = TRUE)
write.csv(filterd_ppi_2 , file = "~/project/input_data.csv", row.names = TRUE)
