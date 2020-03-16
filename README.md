Workflow for Analyzing Single-nuclei-sequencing Pediatric Glioblastoma Cells
====
## Introductions
We analyze the single-nuclei sequencing (sNuc-seq) data of six human brain tumors generated in Dr. Mariella G. Filbin`s lab. We detected astrocytes and neuronal cells by comparing our data to reference datasets. Our primary goals are: finding neurons and astrocytes, identify normal and tumor cells, and check wether any transcripts are secreted. At the current stage, we have found cells similar to neurons and astrocytes, and apply copy number variation analysis to differentiate normal and tumor cells. 
  
## Quantification
We apply HISAT2 (v2.1.0) for alignment, and StringTie (v2.0.3) for transcript assembly and quantification.

To do it, run following R file to generate .sh files :
```{r }
/users/whou/alsf_filbin/quantification/code/01_make_sh.R 
```

And then execute all  the .sh files genereated.

## Data Processing
Gene expression of all cells from all plates of all samples are concatenated as a gene by cell matrix (47127 genes $\times$ 1152 cells).
```{r }
/users/whou/alsf_filbin/data/code/01_logTPM.R
```
We retain cells with at least 1000 genes expressed and alignment rate > 50% (762 cells retained).
We retain genes with log2-scaled TPM > 0.1 in at least 1% of cells (27556 genes retained).

```{r }
/users/whou/alsf_filbin/data/code/02_logTPM_filtered.R
```

Each sample is normalize by library size (by SCRAN).
For each gene, a trend is fitted to the variance against the mean. The fitted value of this trend represents technical variability, for example due to sequencing at a given mean, under the assumption that each cell has the same amount of spike-in RNA and therefore the differences in observed expression account for measurement errors.
SCRAN is applied to decompose the gene-specific variance into biological and technical components by interpolating the fitted trend at the mean log-count for each gene.
Highly variable genes (HVGs) are identified as those with total variance greater than technical variance.
Principal component analysis (PCA) is performed using HVGs.
Top 15 Principal components (PCs) are used to perform further dimension reduction using UMAP for visualization.
```{r }
/users/whou/alsf_filbin/data/code/03_pca_umap.R
```
## Celltypes Identification 
SingleR is applied to identify the celltypes by comparing the cells with two reference datasets: La Manno (2016, Cell) dataset and  Nowakowski (2017, Science) dataset.   
To do it, run the following R file:
```{r }
/users/whou/alsf_filbin/identifyCellType/code/01_identify_celltype.R
```
## Integration and Visualization with Atlas
In this analysis, we identify the cell types of each cell using La Manno (2016, Cell) datasets as the reference.
Since the La Manno dataset consists of mouse and human cells sequenced in fludigm c1 protocol, and Filbin`s data is human cells, we apply Seurat (v3.0) to align these datasets. 
As a data integration method, Seurat utilizes canonical correlation analysis (CCA) to identify correlated gene modules that are present in both datasets, and then mutual nearest neighbors (MNN) of the cells were identified as ``anchors`` to integrate the datasets.
Each sample is integrated with La Manno data separately.

To do it, run the following R file:
```{r }
/users/whou/alsf_filbin/data/integrate_each_to_atlas_rmDEG/code/01_integrate.R
/users/whou/alsf_filbin/data/integrate_each_to_atlas_rmDEG/code/02_umap.R
/users/whou/alsf_filbin/data/integrate_each_to_atlas_rmDEG/code/03_plot_umap.R
```

  
## Copy Number Variation Analysis
We apply inferCNV (v.1.3.4) to compare our data to the normal human cells in La Manno datasets. 
To do it, run the following R file:
```{r }
/users/whou/alsf_filbin/data/infercnv/01_make.R
/users/whou/alsf_filbin/data/infercnv/02_run_infercnv.R
```
  
## Contact the Author
Author: Wenpin Hou

Report bugs and provide suggestions by sending email to:

Maintainer: Wenpin Hou (whou10@jhu.edu)

Or open a new issue on this Github page
