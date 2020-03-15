<<<<<<< HEAD
library(Seurat)
library(data.table)
i = as.numeric(commandArgs(trailingOnly = T)[[1]][1])
source('/dcl02/hongkai/data/whou/resource/function.R')
homolog = fread('/dcl02/hongkai/data/whou/alsf_filbin/homologgene/homologene.data.txt', data.table=F)
homolog = homolog[homolog[,2] %in% c('9606','10090'),]
num <- sapply(unique(homolog[,1]),function(i) {
  sum(homolog[homolog[,1]==i,2]=='9606')==1 & sum(homolog[homolog[,1]==i,2]=='10090')==1
})
tar <- unique(homolog[,1])[which(num)]

t1 <- homolog[homolog[,1] %in% tar & homolog[,2]=='9606',]
t2 <- homolog[homolog[,1] %in% tar & homolog[,2]=='10090',]
identical(t1[,1],t2[,1])
gn <- cbind(t1[,4],t2[,4])
gn[,1] <- toupper(gn[,1])
gn[,2] <- toupper(gn[,2])

countfunc <- function(mat, platform='fluidigm', normalize = TRUE,species='human', log.transform.only = FALSE) {
  row.names(mat) <- toupper(row.names(mat))
  rownames(mat) <- gsub('_LOC.*','',rownames(mat))
  if (normalize){
    if (platform == '10x'){
      mat = log2CPM_from_10x_count(mat)
    } else {
      mat = log2TPM_from_fludigm_count(mat)
    } 
  }
  if (log.transform.only){
    mat = log2(mat + 1)
  }
  mat <- mat[!duplicated(rownames(mat)),]
  if (species=='human') {
    mat <- mat[intersect(row.names(mat),gn[,1]),]
  } else {
    mat <- mat[intersect(row.names(mat),gn[,2]),]
    row.names(mat) <- gn[match(row.names(mat),gn[,2]),1]
  }
  set.seed(12345)
  if (ncol(mat) > 10000) mat <- mat[,sample(1:ncol(mat),10000)]
  rmDupGenesNameOnly(mat, rownames(mat))
}

## load atlas
dlist = list()
print('8')
dlist[['2016_La:ES']] <- countfunc(readRDS('/dcl02/hongkai/data/whou/alsf_filbin/atlases/LaManno_dataset/proc/expr/ES_fulldataset.rds'), species='human')
print('9')
dlist[['2016_La:Human_Embryo']] <- countfunc(readRDS('/dcl02/hongkai/data/whou/alsf_filbin/atlases/LaManno_dataset/proc/expr/Human_Embryo_fulldataset.rds'), species='human')
print('10')
dlist[['2016_La:iPSC']] <- countfunc(readRDS('/dcl02/hongkai/data/whou/alsf_filbin/atlases/LaManno_dataset/proc/expr/iPSC_fulldataset.rds'), species='human')
print('11')
dlist[['2016_La:Mouse_Embryo']] <- countfunc(readRDS('/dcl02/hongkai/data/whou/alsf_filbin/atlases/LaManno_dataset/proc/expr/Mouse_Embryo_fulldataset.rds'), species='mouse')
dlist[['2016_La:MouseAdult']] <- countfunc(readRDS('/dcl02/hongkai/data/whou/alsf_filbin/atlases/LaManno_dataset/proc/expr/MouseAdultDAMoleculeCounts.rds'), species='mouse')

## load data
setwd('/users/whou/')
d = readRDS('./alsf_filbin/data/processed/logTPM_filtered.rds')
rownames(d) = gsub('_.*','',rownames(d))
d = rmDupGenesNameOnly(d, rownames(d))
asp = sub('-.*','',colnames(d))
sp = unique(asp)[i]
rownames(d) = sub('_.*','', rownames(d))
dlist[[sp]] = d[,asp==sp]

gn <- sapply(dlist,row.names)
gn <- table(unlist(gn))
gn <- names(gn)[gn==length(dlist)]
for (i in 1:length(dlist)) dlist[[i]] <- dlist[[i]][gn,]

for (i in names(dlist)) {
  print('findvariablefeatures...')
  print(i)
  dlist[[i]] <- FindVariableFeatures(CreateSeuratObject(dlist[[i]],project=i))
}
print('18')
d.anchors <- FindIntegrationAnchors(object.list = dlist,k.filter=20)
print('19')
d.integrated <- IntegrateData(anchorset = d.anchors)
print('20')
dir.create(paste0('./alsf_filbin/data/integrate_each_to_atlas_rmDEG/result/',sp),showWarnings = F, recursive = T)
saveRDS(d.integrated,paste0('./alsf_filbin/data/integrate_each_to_atlas_rmDEG/result/',sp,'/integrated.rds'))

=======
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

## diffgene
diffgene = readRDS(paste0('./alsf_filbin/data/integrate_each_to_atlas_rmDEG/result/',sp,'/diffgene.rds'))
diffgene = diffgene[!is.na(diffgene[,2]),]
dg = as.character(diffgene[diffgene[,2]>=0.05,'Gene'])
d = d[dg,]
m = m[dg,]
dn = dimnames(m)
m = as.numeric(m)
m = matrix(m,nrow=length(dg))
dimnames(m) = dn

## prepare Seurat object
colnames(m) = paste0('ref;',colnames(m))
colnames(d) = paste0('glioma;',colnames(d))
dlist = list()
dlist[['ref']] = m
dlist[['glioma']] = d[,grep(sp, asp)]
for (i in names(dlist)) {
  dlist[[i]] <- FindVariableFeatures(CreateSeuratObject(dlist[[i]],project=i))
}
reference.vec = c('ref')
set.seed(12345)
d.anchors <- FindIntegrationAnchors(object.list = dlist, reference = grep('ref',names(dlist)), k.filter = 20)
d.integrated <- IntegrateData(anchorset = d.anchors)
dir.create(paste0('./alsf_filbin/data/integrate_each_to_atlas_rmDEG/result/',sp),showWarnings = F, recursive = T)
saveRDS(d.integrated, paste0('./alsf_filbin/data/integrate_each_to_atlas_rmDEG/result/',sp,'/integrated.rds'))
>>>>>>> 493eb965c5992a7a722bc79e959aa16119876f2b

