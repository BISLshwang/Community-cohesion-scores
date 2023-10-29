# Community cohesion looseness in genetic networks reveals individualized drug targets and resistance for precision oncology
Community cohesion originally introduced in social science represents how closely people in a particular community interact with each other. It has been known that the social interactions reduce risk of diseases and mortality. Herein, we propose that the concept of community cohesion could be applied between genes as well for precision oncology to overcome limitations of single gene-based biomarkers. We develop community cohesion scores which capture unique characteristics of the community cohesion looseness in each individualized genetic network. It allows precise quantification of the community ability to retain the interactions between the genes and their cellular functions. Using breast cancer as proof-of-concept study to illustrate our method, we show that the community cohesion scores can be used as superior therapeutic and prognostic biomarkers to predict individualized drug targets and resistance than the previous approaches. Our method enables to open a new horizon for the biomarker development in precision oncology.

## 1. Dataset
All data (gene expression profiles) used in the paper are publicly available. 
For example, 
1) GTEx (GTEx Consortium. "The GTEx Consortium atlas of genetic regulatory effects across human tissues." Science 369.6509 (2020): 1318-1330.)
  - https://gtexportal.org/home/
2) TCGA (Weinstein, John N., et al. "The cancer genome atlas pan-cancer analysis project." Nature genetics 45.10 (2013): 1113-1120.)
  - https://portal.gdc.cancer.gov/
3) GEO (Clough, Emily, and Tanya Barrett. "The gene expression omnibus database." Statistical Genomics: Methods and Protocols (2016): 93-110.) 
  - https://www.ncbi.nlm.nih.gov/geo/
4) GDSC (Yang, Wanjuan, et al. "Genomics of Drug Sensitivity in Cancer (GDSC): a resource for therapeutic biomarker discovery in cancer cells." Nucleic acids research 41.D1 (2012): D955-D961.)
  - https://www.cancerrxgene.org/


## 2. Source code
1) Code for constructing normal tissue genetic networks
   - WCGNA.R
2) Code for estimating individaulized genetic networks and quantifying community cohesion scores
   - community_cohesion_scores.py

Instruction to use the source code:
1) Download the code.
2) Make sure you are in the main directory where the code is.
3) Run the following.

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
