# Community cohesion looseness in genetic networks reveals individualized drug targets and resistance 
Community cohesion originally introduced in social science represents how closely people in a particular community interact with each other. It has been known that the social ties reduce risk of diseases and mortality. Here, we propose that the concept of community cohesion could be applied between genes as well to discover therapeutic and prognostic biomarkers for precision oncology. We develop community cohesion scores which capture unique characteristics of the community cohesion looseness in each individualized genetic network. It allows precise quantification of the community ability to retain the normal interactions between the genes and their cellular functions. Using breast cancer as proof-of-concept study to illustrate our method, we show that the community cohesion scores are superior to single-gene based biomarkers in terms of identifying individualized drug targets and resistance. Our method enables to open a new horizon for the biomarker development in precision oncology.

## 1. Dataset
All data (gene expression profiles) used in the paper are publicly available. 
For example, 
1) GTEx (GTEx Consortium. "The GTEx Consortium atlas of genetic regulatory effects across human tissues." Science 369.6509 (2020): 1318-1330, https://gtexportal.org/home/)
2) TCGA (Weinstein, John N., et al. "The cancer genome atlas pan-cancer analysis project." Nature genetics 45.10 (2013): 1113-1120., https://portal.gdc.cancer.gov/)
3) GEO (Clough, Emily, and Tanya Barrett. "The gene expression omnibus database." Statistical Genomics: Methods and Protocols (2016): 93-110., https://www.ncbi.nlm.nih.gov/geo/)
4) GDSC (Yang, Wanjuan, et al. "Genomics of Drug Sensitivity in Cancer (GDSC): a resource for therapeutic biomarker discovery in cancer cells." Nucleic acids research 41.D1 (2012): D955-D961., https://www.cancerrxgene.org/)


## 2. Source code
1) Code for constructing the normal tissue genetic networks
   - WGCNA.R
2) Code for estimating the individaulized genetic networks and quantifying the community cohesion scores
   - community_cohesion_scores.py

Instruction to use the source code:
1) Download the code.
2) Make sure you are in the main directory where the code is.
3) Run the following.
   - Prepare 4 types of input files (e.g., we provide toy examples files in /examples)
      - df_comm.csv: file providing the genes and the corresponding genetic community
      - df_weight.csv: file providing the network edges and their weights (we only provides the weights of one community due to the uploded file size limitation)
      - control_samples_breast.csv: file providing the gene expression profiles of control samples
      - case_samples_breast_cancer.csv: file providing the gene expression profiles of case samples (at least one sample)

```
python3 community_cohesion_scores.py examples/df_comm.csv examples/df_weight.csv examples/control_samples_breast.csv examples/case_samples_breast_cancer.csv
```

## 3. Libaray
1) python==3.7.13
2) numpy==1.21.5
3) pandas==1.3.5
4) scipy==1.7.3
5) networkx==2.5.1
6) tqdm==4.62.1
