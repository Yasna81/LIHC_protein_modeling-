# LIHC_protein_modeling-
  ###  🧬  Protein Interaction Prediction in Hepatocellular Carcinoma using Transcriptomics and Functional Similarity.


  
  ----- 
  ####  🚀 Motivation :


This project was based on the idea that transcriptomic data can indirectly reflect functional protein associations, especially when combined with biological context such as Gene Ontology (GO) similarity and prior knowledge from STRINGdb. By leveraging gene co-expression and semantic similarity, we aimed to build a predictive model capable of identifying potential protein-protein interactions (PPIs) specific to liver cancer biology. This model could help narrow down candidate interactions for wet-lab validation.


  
  ----
  ####  💻 Methods:

We first downloaded RNA-Seq data from TCGA-LIHC (hepatocellular carcinoma) and identified differentially expressed genes (DEGs) using DESeq2. We applied strict thresholds (padj < 0.05, |log2FoldChange| > 1) for significant genes. Next, we extracted protein-protein interaction data from STRINGdb (v11.5, Homo sapiens) and filtered it down to high-confidence interactions. Starting with over 517,146 pairs (combined score > 700), we tightened the criteria to just 3,314 pairs (combined score > 900) to focus on the strongest signals.

To assign biologically meaningful labels, we calculated Gene Ontology (GO) similarity (Biological Process ontology) for all pairs using the GOSemSim package. Pairs with go_sim ≥ 0.6 were labeled as interacting (1), while those with go_sim ≤ 0.4 were labeled as non-interacting (0). Pairs with intermediate similarity (0.4 < go_sim < 0.6) were excluded to minimize ambiguity, leaving 1,378 pairs for further analysis.

For each PPI, we computed co-expression Pearson correlation of RNA-seq expression levels, which was calculated only on training data to prevent data leakage. We trained a logistic regression model (GLM) to predict interactions. Co-expression correlation and PPI combined score were used as predictors, while the response was binary labels ((1) = 519 / (0) = 859) from GO semantics.

The final model was evaluated by both biological validation and an independent test set (30% held-out data). Gene Set Enrichment Analysis (GSEA) of top-predicted genes was conducted with the clusterProfiler package in R.
____
####  🏁 Results :

From the LIHC RNA-seq data, we identified 1,232 significant differentially expressed genes using DESeq2; these genes served as the candidate set for potential PPIs. A heatmap of the top 50 significant genes is illustrated below.


<img width="1181" height="682" alt="HEATMAP" src="https://github.com/user-attachments/assets/3b4370c0-28d1-4ae1-8f09-b42bffe21ed8" />


Figure 1



The top upregulated gene included COX7B2 (log2FC = 3.9), while BMP10 (log2FC = -3.6) was among the most downregulated.


<img width="555" height="320" alt="vocano_enh" src="https://github.com/user-attachments/assets/b59ec3b6-456f-47b3-858a-f3b1e0edd4f4" />


Figure 2

Using STRINGdb high-confidence interactions (combined_score > 900) and GO semantic similarity, we labeled 519 pairs as interacting (go_sim > 0.6) and 859 as non-interacting (go_sim < 0.4). The logistic regression model trained on co-expression (Pearson) and STRING combined scores achieved an AUC of 0.77 on the test set (Figure 3), with a sensitivity of 0.82, specificity of 0.6, and accuracy of 0.72. The top predicted novel interaction was RPL38 with RPL18, with a predicted probability of 0.82 and a combined score of 999.


<img width="555" height="320" alt="ROC curve" src="https://github.com/user-attachments/assets/81ded011-1a83-48c1-8321-745118aa7a24" />


Figure 3

To assess the biological validation of the model, we conducted gene set enrichment analysis of genes from top-predicted PPIs. KEGG enrichment analysis (Figure 4) suggests cancer-related pathways, including cell cycle regulation, spliceosome function, ubiquitin-mediated proteolysis, and TGF-β signaling. Additionally, Coronavirus disease - COVID-19 and neurodegenerative pathways were identified, which share similar pathways with cancers.

<img width="700" height="563" alt="kegg" src="https://github.com/user-attachments/assets/99667b6d-082b-4814-a565-0df8996b5140" />


Figure 4

Gene Ontology (GO) enrichment analysis suggests transcriptional and post-translational regulation, RNA processing, and protein synthesis - processes that are frequently dysregulated in cancer. These findings suggest that our interaction-based classifier captures biologically meaningful and cancer-associated features.


<img width="700" height="563" alt="go" src="https://github.com/user-attachments/assets/965c8afb-92cc-4456-b1d5-2cb487ae2f8f" />



Figure 5

The PPI network shown below visualizes the top 20 predicted interacting gene pairs based on our machine learning model. Each node represents a gene (protein), and each edge indicates a predicted interaction with high confidence.




<img width="603" height="494" alt="20-top" src="https://github.com/user-attachments/assets/525e84ac-05e2-4a5e-81b2-28fea79566ee" />





Figure 6

To evalute the biological relevence of the predicted protein-protein interactions we preformed validation usining the STRING database (v11.5), gene pairs predicted by the model were matched against SRING interactions with a combined confidencce score > 900. Out of 351 predicted gene pairs (score >900) , 178  were found in STRING high confidence interaction network, resulting in a precision of **25%**. While this may represent false positive arising from computational limitations ( expression without physical interaction) This could reflect novel specific interactions which needs experimental validation to confirm the biological relevence. 



<img width="526" height="375" alt="evolution" src="https://github.com/user-attachments/assets/c936f7a5-e5b7-47ab-896a-d9533844cf39" />



Figure 7

____

#### 📦 packages used :

 DESeq2 ,STRINGdb ,biomaRt ,GOSemSim, org.Hs.eg.db,pROC,caret,clusterProfiler,igraph,pheatmap,ggplot2,ggrepel,visNetwork


 ▶️ Input and output data are available in .csv

-----

#### 🥇 future enhancments:

