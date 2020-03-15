<<<<<<< HEAD
mat = readRDS('/users/whou/alsf_filbin/data/processed/logTPM_filtered.rds') ##  [1:27556, 1:762] 
data <- mat
suppressMessages(library(scran))
sep <- '-'
a = sapply(colnames(data), function(i) strsplit(i,'-')[[1]][1])
b = sapply(colnames(data), function(i) strsplit(i,'-')[[1]][2])
batch = paste0(a,'-',b)
=======
mat = readRDS('/users/whou/alsf_filbin/data/processed/logTPM_filtered.rds') ## [1:29229, 1:299] 
data <- mat
suppressMessages(library(scran))
sep <- '-'
batch <- sub(paste0(sep,'.*'),'',colnames(data))
>>>>>>> 493eb965c5992a7a722bc79e959aa16119876f2b
id <- order(batch)
data <- data[,id]
batch <- batch[id]
## highly variable genes
cn <- colnames(data)
fit <- trendVar(data)
decomp <- decomposeVar(data,fit)
gs <- row.names(decomp)[decomp[,'total'] > decomp[,'tech']]
<<<<<<< HEAD
data <- data[gs,]  
=======
data <- data[gs,]  ###  [1:13662, 1:299]
>>>>>>> 493eb965c5992a7a722bc79e959aa16119876f2b
scmd <- sapply(1:length(unique(batch)),function(sp) {
  paste0("data[gs,batch==unique(batch)[",sp,"]]")
})

pr = prcomp(t(data),scale. = T)

<<<<<<< HEAD
dir.create('/users/whou/alsf_filbin/data/plot')
=======
dir.create('/home-4/whou10@jhu.edu/scratch/Wenpin/filbin/plot')
>>>>>>> 493eb965c5992a7a722bc79e959aa16119876f2b
pdf('/users/whou/alsf_filbin/data/plot/pca_statistic.pdf',width=12,height=3)
par(mfrow=c(1,4))
results = decomp
plot(results$mean, results$total,xlab='Mean normalized log-expression per gene',ylab='Variance of the normalized log-expression per gene')
o <- order(results$mean)
lines(results$mean[o], results$tech[o], col="red", lwd=2)
plot(results$mean, results$bio,xlab='Mean normalized log-expression per gene',ylab='Biological component of the variance')
plot(results$total, results$tech,xlab='Variance of the normalized log-expression per gene',ylab='Technical component of the variance')
abline(a=0,b=1,col='red')
plot(pr$sdev[1:100],xlab='number of PC', ylab='sd',main='',pch=20)
abline(v=15,col='red')
dev.off()

<<<<<<< HEAD
saveRDS(pr,file='/users/whou/alsf_filbin/data/processed/pr_all.rds')
=======
>>>>>>> 493eb965c5992a7a722bc79e959aa16119876f2b
pr <- pr$x[,1:15]
saveRDS(pr,file='/users/whou/alsf_filbin/data/processed/pr.rds')

## pca
tab = table(batch)
batch2 = paste0(batch,'(',tab[match(batch, names(tab))],')')
library(ggplot2)
<<<<<<< HEAD
library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
pdf('/users/whou/alsf_filbin/data/plot/PCA.pdf',width=6,height=3.5)
ggplot() + geom_point(data=data.frame(pr[,1:2],batch=as.factor(batch2)),aes(x=PC1,y=PC2,col=batch),alpha=0.8, size=0.5) + theme_classic() + xlab('PC1') + ylab('PC2') + theme(legend.title = element_blank())
=======
pdf('/users/whou/alsf_filbin/data/plot/PCA.pdf',width=6,height=4.5)
ggplot() + geom_point(data=data.frame(pr[,1:2],batch=as.factor(batch2)),aes(x=PC1,y=PC2,col=batch),alpha=0.8) + theme_classic() + xlab('PC1') + ylab('PC2')
>>>>>>> 493eb965c5992a7a722bc79e959aa16119876f2b
dev.off()

## umap
library(umap)
set.seed(12345)
u <- umap(pr)$layout
<<<<<<< HEAD
saveRDS(u, file='/users/whou/alsf_filbin/data/processed/umap.rds')
pdf('/users/whou/alsf_filbin/data/plot/umap.pdf',width=6,height=3.5)
ggplot() + geom_point(data=data.frame(umap1=u[,1],umap2=u[,2],batch=as.factor(batch2)),aes(x=umap1,y=umap2,col=batch),alpha=0.8, size=0.5) + theme_classic() + xlab('UMAP1') + ylab('UMAP2')
dev.off()

pdf('/users/whou/alsf_filbin/data/plot/umap_batch.pdf',width=8,height=4)
ggplot() + geom_point(data=data.frame(umap1=u[,1],umap2=u[,2],batch=as.factor(batch2)),aes(x=umap1,y=umap2,col=batch),alpha=0.8, size=0.5) + theme_classic() + facet_wrap(~batch) + theme(legend.position = 'none') + xlab('UMAP1') + ylab('UMAP2')
dev.off()

af <- list.files('/dcl02/hongkai/data/whou/alsf_filbin/align/')
qua <- t(sapply(af,function(f) {
  d <- readLines(paste0('/dcl02/hongkai/data/whou/alsf_filbin/align/',f,'/summary.txt'))
  c(as.numeric(sub(' .*','',d[1])),as.numeric(sub('%.*','',d[length(d)])))
}))
colnames(qua) = c('num.read','alignment_rate')
qua = qua[colnames(mat), ]
library(gridExtra)
pdf('/users/whou/alsf_filbin/data/plot/umap_quality.pdf',width=11,height=2)
p1 <- ggplot() + geom_point(data=data.frame(umap1=u[,1],umap2=u[,2],batch=as.factor(batch2), log10_num.read = log10(qua[,1]), alignment.rate = qua[,2]),aes(x=umap1,y=umap2,col=log10_num.read,),alpha=0.8, size=0.5) + theme_classic() + xlab('UMAP1') + ylab('UMAP2')
p2 <- ggplot() + geom_point(data=data.frame(umap1=u[,1],umap2=u[,2],batch=as.factor(batch2), num.read = qua[,1], alignment.rate = qua[,2]),aes(x=umap1,y=umap2,col=alignment.rate),alpha=0.8, size=0.5) + theme_classic()+ xlab('UMAP1') + ylab('UMAP2')
p3 <- ggplot() + geom_point(data=data.frame(umap1=u[,1],umap2=u[,2],batch=as.factor(batch2), num.read = qua[,1], alignment.rate = qua[,2], sum.expression=colSums(mat)),aes(x=umap1,y=umap2,col=sum.expression),alpha=0.8, size=0.5) + theme_classic()+ xlab('UMAP1') + ylab('UMAP2')
grid.arrange(p1,p2,p3,nrow=1)
dev.off()

pdf('/users/whou/alsf_filbin/data/plot/umap_quality_gene.pdf',width=11,height=2)
p1 <- ggplot() + geom_point(data=data.frame(umap1=u[,1],umap2=u[,2],batch=as.factor(batch2), MBP = mat[grep('^MBP_',rownames(mat)),]),aes(x=umap1,y=umap2,col=MBP),alpha=0.8, size=0.5) + theme_classic()+ xlab('UMAP1') + ylab('UMAP2')
p2 <- ggplot() + geom_point(data=data.frame(umap1=u[,1],umap2=u[,2],batch=as.factor(batch2), PLP1 = mat[grep('^PLP1_',rownames(mat)),]),aes(x=umap1,y=umap2,col=PLP1),alpha=0.8, size=0.5) + theme_classic()+ xlab('UMAP1') + ylab('UMAP2')
p3 <- ggplot() + geom_point(data=data.frame(umap1=u[,1],umap2=u[,2],batch=as.factor(batch2), TF = mat[grep('^TF_',rownames(mat)),]),aes(x=umap1,y=umap2,col=TF),alpha=0.8, size=0.5) + theme_classic()+ xlab('UMAP1') + ylab('UMAP2')
grid.arrange(p1,p2,p3,nrow=1)
dev.off()

m1 = mat[,u[,2]>5]
m2 = mat[,u[,2]<5]
fc = rowMeans(m1) - rowMeans(m2)
sort(fc,decreasing = T)[1:20]
summary(rowMeans(m1[grepl('^MT-',rownames(m1)), ]))
summary(rowMeans(m2[grepl('^MT-',rownames(m2)), ]))
=======
pdf('/users/whou/alsf_filbin/data/plot/umap.pdf',width=6,height=4.5)
ggplot() + geom_point(data=data.frame(umap1=u[,1],umap2=u[,2],batch=as.factor(batch2)),aes(x=umap1,y=umap2,col=batch),alpha=0.8) + theme_classic()
dev.off()
>>>>>>> 493eb965c5992a7a722bc79e959aa16119876f2b
