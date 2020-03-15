i = as.numeric(commandArgs(trailingOnly = T)[[1]])
library(Seurat)
## load data
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/')
d = readRDS('./alsf_filbin/data/processed/logTPM_filtered.rds')
asp = sub('-.*','',colnames(d))
sp = unique(asp)[i]
rownames(d) = sub('_.*','', rownames(d))

## load atlas
m = read.table('./alsf_filbin/data/Kriegstein_Human_cortex/exprMatrix.tsv',as.is=T)
m = as.matrix(m)
rownames(m) = m[,1]
m = m[,-1]
colnames(m) = m[1,]
m = m[-1,]

## intersect genes
intergene = intersect(rownames(m), rownames(d))
d = d[intergene,]
m = m[intergene,]
dn = dimnames(m)
m = as.numeric(m)
m = matrix(m,nrow=length(intergene))
dimnames(m) = dn
d_s = d[,grep(sp, asp)]
## prepare Seurat object
pval <- sapply(1:nrow(m), function(g){
  wilcox.test(m[g,], d_s[g,])$p.value
})
fdrv <- p.adjust(pval,method='fdr')
df = data.frame(Gene=rownames(m), fdr = fdrv)
dir.create(paste0('./alsf_filbin/data/integrate/result/',sp),showWarnings = F, recursive = T)
saveRDS(df, paste0('./alsf_filbin/data/integrate/result/',sp,'/diffgene.rds'))
