[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6423597.svg)](https://doi.org/10.5281/zenodo.6423597)

﻿﻿﻿﻿﻿﻿﻿﻿﻿﻿﻿﻿﻿Introduction
------------------------

README file for codes to reproduce the three downstream analyses (such as *differential expression analysis*, *cell clustering analysis* and *pseudotime analysis*) in the papers **"scTSSR: gene expression recovery for single-cell RNA
sequencing using two-side sparse self-representation"** and **"scTSSR2: imputing dropout events for single-cell RNA sequencing using fast two-side
self-representation"** . The methods `scTSSR` and `scTSSR2` developed in the papers can be found at: [https://github.com/Zhangxf-ccnu/scTSSR](https://github.com/Zhangxf-ccnu/scTSSR); [https://github.com/Zhangxf-ccnu/scTSSR2](https://github.com/Zhangxf-ccnu/scTSSR2).
To reproduce the *comparison with smRNA FISH data* and *down-sampling* experiments, please refer to Huang et al. (2018): [SAVER-paper](https://github.com/mohuangx/SAVER-paper).


Contents of this archive
------------------------
This archive contains:   
 
(1) **Datasets**: subdirectory that contains four preprocessed datasets: bulk data and single-cell data of **H1\_DEC** (such as **H1\_DEC\_bulk** and **H1\_DEC\_sc**), **sc\_celseq2\_5cl\_p1** and **Deng** datasets, which can be used to reproduce the three downstream analyses (such as *differential expression analysis*, *cell clustering analysis* and *pseudotime analysis*) respectively.   
Note that **H1\_DEC** data has been preprocessed using Seurat v3.2 to contain 2,000 highly variable genes, and **sc\_celseq2\_5cl\_p1** data has been preprocessed following the method of [Hou et al.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02132-x), and **Deng** dataset has been preprocessed by filtering out genes expressed in less than 10% of cells.


(2) **DE\_analysis.R**: Using `edgeR` to identify differential expression genes (DEGs) between pairs of cell subpopulations.  And `AUC` and `Spearman correlation coefficient` are computed as the criteria.

(3) **Cell\_clustering.R**: Using `SC3` and `Seurat` to carry out cell clustering analysis. And `ARI` is used to evaluate the consistency between the results of `SC3` or `Seurat` and the reference labels of cells. 

(4) **Trajectory\_inference.R**: Using `TSCAN` and `Monocle2` to carry out pseudotime analysis. The function returns `POS` and `Kendall's rank correlation` scores.  



Contact
------------------------
Please do not hesitate to contact Miss **Ke Jin** [kej13@mails.ccnu.edu.cn](kej13@mails.ccnu.edu.cn) or Dr. **Xiao-Fei Zhang** [zhangxf@mail.ccnu.edu.cn](zhangxf@mail.ccnu.edu.cn) to seek any clarifications regarding any contents or operation of the archive.















































