# LIHC_protein_modeling-
  ###  🧬  Protein Interaction Prediction in Hepatocellular Carcinoma using Transcriptomics and Functional Similarity.


  
  ----- 
  ####  🚀 Motivation :


This project was based on the idea that transcriptomic data can indirectly reflect functional protein associations, especially when combined with biological context such as Gene Ontology (GO) similarity and prior knowledge from STRINGdb. By leveraging gene co-expression and semantic similarity, we aimed to build a predictive model capable of identifying potential protein-protein interactions (PPIs) specific to liver cancer biology. This model could help narrow down candidate interactions for wet-lab validation.


  
  ----
  ####  💻 Methods:

we first downloaded RNA-Seq data from TCGA-LIHC (hepatocellular carcinoma) and identified differentially expressed genes (DEGs) using DESeq2. We applied strict thresholds (padj < 0.05, |log2FoldChange| > 1) for significant genes. Next, we extracted protein-protein interaction data from STRINGdb (v11.5, Homo sapiens) and filtered it down to high-confidence interactions. Starting with over 517146  pairs (combined score > 700), we tightened the criteria to just 3,314 pairs (combined score > 900) to focus on the strongest signals.

To assign biologically meaningful labels, we calculated Gene Ontology (GO)similarity (Biological Process ontology) for all pairs using GOSemSim package. Pairs with go_sim >= 0.6 were labeld as interacting (1) while those with go_sim <= 0.4 were label as non-interacting (0), pairs with intermediate similarity (0.4 < go_sim < 0.6) were excluded to minimize ambiguity leaving 1378 pairs for further analysis.

For each PPI we computed CO-expression  Pearson correlation of RNA-seq expression levels ,which was calculated only on training data to prevent data leakage. We trained a logistic regression model (glm) to predict interactions. Co-expression correlation and  PPI combined-score  were used as predictors while response was binary lables ((1) = 519 /(0) = 859) from GO semantics.

Finall Model was evaluted by both biological validation and Independent test set (30% held-out data).Gene Set enrichmet analysis (GSEA) of top-predicted genes was conducted by clusterProfiler package in R 

____
####  🏁 Results :

From the LIHC RNA-seq data, we identified 1232 significant differentially expressed genes using Desq2, these genes served as the candidate set for potential PPIs. A heatmap of top 50 significant genes are illustreted below.
<img width="1181" height="682" alt="HEATMAP" src="https://github.com/user-attachments/assets/3b4370c0-28d1-4ae1-8f09-b42bffe21ed8" />




Top upregulated gene included COX7B2 (log2FC = 3.9 ) while BMP10 (log2FC = -3.6 ) was among the most downregulated.
<img width="555" height="320" alt="vocano_enh" src="https://github.com/user-attachments/assets/b59ec3b6-456f-47b3-858a-f3b1e0edd4f4" />


Figure_2

Using STRINGdb high_confidence interactions (combined_score > 900) and GO semantic similarity, we labeled 519 pairs as interacting (go_sim > 0.6) and 859 as non-interacting (go_sim < 0.4). Logistic regression model trained on co-expresion (Pearson) and STRING combined scores achieved AUC of 0.77 on test set(figureB), Sensitivity of 0.82 , Specificity of 0.6 and Accurecy of 0.72. Top predicted novel interaction was RPL38 and RPL18 with predicted-probability of 0.82 and combined-score of 999 .


<img width="555" height="320" alt="ROC curve" src="https://github.com/user-attachments/assets/81ded011-1a83-48c1-8321-745118aa7a24" />


Figure_3

To assess the Biological Validation of the model we conducted gene set enrichment analysis of genes from top-predicted PPIs. KEGG enrichment analysis(figurec) suggests cancer related pathwyas including cell cycle regulation, splicesome function, ubiquitin-mediated proteolysis and TGF-B signaling. Also Coronavirus disease-Covid19 and Neurodegenrative pathways that share similiar pathways with cancrs. 

<img width="700" height="563" alt="kegg" src="https://github.com/user-attachments/assets/99667b6d-082b-4814-a565-0df8996b5140" />


Figure_4






Gene Ontology (GO) enrichment analysis suggests transcriptional and post-translational regulation, RNA processing and protien synthesis that are frequently dysregulated in cancer. These findings suggests that our interaction-based classifier captures biologically meaningful and cancer-associated features.


<img width="700" height="563" alt="go" src="https://github.com/user-attachments/assets/965c8afb-92cc-4456-b1d5-2cb487ae2f8f" />

Figure_5

                                                          
The PPI network shown below visulizes the 20 top predicted interacting gene pairs based on our machine learning model. Each node represents a gene(protien) and each edge indicates a predicted interaction with high confidence.




<img width="603" height="494" alt="20-top" src="https://github.com/user-attachments/assets/525e84ac-05e2-4a5e-81b2-28fea79566ee" />





Figure_6
____

#### 📦 packages used :


-----

#### 🥇 future enhancments:
