library(data.table)
gtf <- fread('/dcl02/hongkai/data/whou/resource/gtf/grch38.gtf',data.table=F)
gtf <- gtf[gtf[,3]=='gene',]
gn <- sub('\";.*','',sub('.*; gene_name \"','',gtf[,9]))
d <- readRDS('/users/whou/alsf_filbin/data/processed/logTPM_filtered.rds')
rownames(d) = sub('_.*','',rownames(d))
d = d[!duplicated(rownames(d)), ]
source('/users/whou/alsf_filbin/data/integrate_each_to_atlas_rmDEG/code/function.R')
a1 <- countfunc(readRDS('/dcl02/hongkai/data/whou/alsf_filbin/atlases/LaManno_dataset/proc/expr/ES_fulldataset.rds'), species='human')

a2 <- countfunc(readRDS('/dcl02/hongkai/data/whou/alsf_filbin/atlases/LaManno_dataset/proc/expr/Human_Embryo_fulldataset.rds'), species='human')

a3 <- countfunc(readRDS('/dcl02/hongkai/data/whou/alsf_filbin/atlases/LaManno_dataset/proc/expr/iPSC_fulldataset.rds'), species='human')


intgene = intersect(rownames(a1), intersect(rownames(a2), rownames(a3)))
a = cbind(a1[intgene,], a2[intgene,], a3[intgene,])
# set.seed(12345)
# a <- a[,sample(colnames(a),500)]
int <- intersect(row.names(d),row.names(a))
tmp <- cbind(d[int,],a[int,])
tmp <- tmp[row.names(tmp) %in% gn,]
write.table(tmp,file='/users/whou/alsf_filbin/infercnv/genematrix.txt',quote=F,sep='\t')

# meta <- data.frame(cell = colnames(tmp),study = ifelse(grepl('2016_La',colnames(tmp)),'Normal','FilbinData'),stringsAsFactors = F)
meta <- data.frame(cell = colnames(tmp), sample = ifelse(grepl('2016_La',colnames(tmp)),'Normal',sub('-.*','',colnames(tmp))),stringsAsFactors = F)

write.table(meta,file='/users/whou/alsf_filbin/infercnv/meta.txt',quote=F,sep='\t',col.names=F,row.names=F)

gr <- data.frame(gene = row.names(tmp),gtf[match(row.names(tmp),gn),c(1,4,5)],stringsAsFactors = F)
colnames(gr) = c('gene','chr','start','end')
write.table(gr,file='/users/whou/alsf_filbin/infercnv/gr.txt',quote=F,sep='\t',col.names=F,row.names=F)


