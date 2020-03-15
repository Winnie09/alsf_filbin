i = as.numeric(commandArgs(trailingOnly = T)[[1]])
library(Seurat)
<<<<<<< HEAD
setwd('/users/whou/')
af = list.files('./alsf_filbin/data/integrate_each_to_atlas_rmDEG/result/')
u = readRDS(paste0('./alsf_filbin/data/integrate_each_to_atlas_rmDEG/result/',af[i],'/umap/umap.rds'))
cn = rownames(u)
type = sub(':.*','',rownames(u))
type = sub('-.*','',type)

a1 = readRDS('/dcl02/hongkai/data/whou/alsf_filbin/atlases/LaManno_dataset/proc/meta/ES_fulldataset.rds')
a2 = readRDS('/dcl02/hongkai/data/whou/alsf_filbin/atlases/LaManno_dataset/proc/meta/Human_Embryo_fulldataset.rds')
a3 = readRDS('/dcl02/hongkai/data/whou/alsf_filbin/atlases/LaManno_dataset/proc/meta/iPSC_fulldataset.rds')
a4 = readRDS('/dcl02/hongkai/data/whou/alsf_filbin/atlases/LaManno_dataset/proc/meta/MouseAdultDAMoleculeCounts.rds')
# a5 = readRDS('/dcl02/hongkai/data/whou/alsf_filbin/atlases/LaManno_dataset/proc/meta/Mouse_Embryo_fulldataset.rds')

meta = rbind(a1,a2,a3,a4)
ct = c(meta[match(cn[type=='2016_La'], meta$cell),'celltype' ], rep(af[i],sum(type!='2016_La')))
region = c(meta[match(cn[type=='2016_La'], meta$cell),'location' ], rep(af[i],sum(type!='2016_La')))
=======
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/')
af = list.files('./alsf_filbin/data/integrate_each_to_atlas_rmDEG/result/')
u = readRDS(paste0('./alsf_filbin/data/integrate_each_to_atlas_rmDEG/result/',af[i],'/umap/umap.rds'))
cn = sub('.*;','',rownames(u))
type  = sub(';.*','',rownames(u))
meta = read.csv('./alsf_filbin/data/Kriegstein_Human_cortex/meta.csv',as.is=T)
ct = c(meta[match(cn[type=='ref'], meta$Cell),'WGCNAcluster' ], rep(af[i],sum(type!='ref')))
ct[grepl('EN-V1', ct)] <- 'EN-V1'
ct[grepl('IN-CTX-CGE',ct)] <- 'IN-CTX-CGE'
ct[grepl('IN-CTX-MGE',ct)] <- 'IN-CTX-MGE'
ct[grepl('IPC-div',ct)] <- 'IPC-div'
ct[grepl('IPC-nEN',ct)] <- 'IPC-nEN'
ct[grepl('MGE-IPC',ct)] <- 'MGE-IPC'
ct[grepl('MGE-RG',ct)] <- 'MGE-RG'
ct[grepl('nEN-early',ct)] <- 'nEN-early'
ct[grepl('RG-div',ct)] <- 'RG-div'
ct[grepl('EN-PFC',ct)] <- 'EN-PFC'
ct[grepl('nIN',ct)] <- 'nIN'
ct[grepl('U1',ct)] <- 'U'
ct[grepl('U2',ct)] <- 'U'
ct[grepl('U3',ct)] <- 'U'
ct[grepl('U4',ct)] <- 'U'
region = c(meta[match(cn[type=='ref'], meta$Cell),'RegionName' ], rep(af[i],sum(type!='ref')))
>>>>>>> 493eb965c5992a7a722bc79e959aa16119876f2b
pd = data.frame(cell = cn, type = type, celltype=ct, region = region, umap1 = u[,1], umap2=u[,2])
pd = pd[sapply(as.character(pd$celltype), nchar)!=0, ]
pd = pd[!is.na(as.character(pd$celltype)), ]
pd$celltype = factor(as.character(pd$celltype),levels=unique(pd$celltype))
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(gridExtra)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
ct_col = getPalette(length(unique(pd$celltype))+1)[1:length(unique(pd$celltype))]
names(ct_col) = unique(pd$celltype)
region_col = getPalette(length(unique(pd$region))+1)[1:length(unique(pd$region))]
names(region_col) = unique(pd$region)
ct_col[names(ct_col)==af[i]] <-  region_col[names(region_col)==af[i]] <-'cyan'

<<<<<<< HEAD
=======
## mark GBM 'celltype' as 'GBM'
>>>>>>> 493eb965c5992a7a722bc79e959aa16119876f2b
plotdir = './alsf_filbin/data/integrate_each_to_atlas_rmDEG/plot/'
dir.create(plotdir, showWarnings = F, recursive = T)
pdf(paste0(plotdir,af[i],'.pdf'),height=9,width=13)
ggplot() + geom_point(data=pd,aes(x=umap1,y=umap2,col=celltype),alpha=0.8,size=0.5) + 
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'right',legend.title = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
  scale_color_manual(values=ct_col)+
  facet_wrap(~celltype) + theme(legend.position = 'none')
dev.off()
<<<<<<< HEAD

pdf(paste0(plotdir,af[i],'_one.pdf'),height=4,width=4)
ggplot() + geom_point(data=pd,aes(x=umap1,y=umap2,col=celltype),alpha=0.8,size=0.5) + 
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'right',legend.title = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
  scale_color_manual(values=ct_col)+
  theme(legend.position = 'none')
dev.off()

=======
>>>>>>> 493eb965c5992a7a722bc79e959aa16119876f2b
