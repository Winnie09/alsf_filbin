library(data.table)
gtf <- fread('/home-4/zji4@jhu.edu/scratch/resource/gtf/grch38.gtf',data.table=F)
gtf <- gtf[gtf[,3]=='gene',]
gn <- sub('\";.*','',sub('.*; gene_name \"','',gtf[,9]))
d <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/alsf_filbin/data/processed/logTPM_filtered.rds')
a <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Fan_CellResearch/proc/expr.rds')
a <- a[,sample(colnames(a),500)]
a <- log2(a + 1)
row.names(d) <- sub('_.*','',row.names(d))
int <- intersect(row.names(d),row.names(a))
tmp <- cbind(d[int,],a[int,])
tmp <- tmp[row.names(tmp) %in% gn,]
write.table(tmp,file='/home-4/whou10@jhu.edu/scratch/Wenpin/alsf_filbin/data/infercnv/genematrix.txt',quote=F,sep='\t')

meta <- data.frame(colnames(tmp),ifelse(grepl('2018_Fan_CellResearch',colnames(tmp)),'Normal','filbin'),stringsAsFactors = F)
write.table(meta,file='/home-4/whou10@jhu.edu/scratch/Wenpin/alsf_filbin/data/infercnv/meta.txt',quote=F,sep='\t',col.names=F,row.names=F)

gr <- data.frame(row.names(tmp),gtf[match(row.names(tmp),gn),c(1,4,5)],stringsAsFactors = F)
write.table(gr,file='/home-4/whou10@jhu.edu/scratch/Wenpin/alsf_filbin/data/infercnv/gr.txt',quote=F,sep='\t',col.names=F,row.names=F)

