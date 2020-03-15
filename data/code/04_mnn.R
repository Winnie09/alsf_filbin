mat = readRDS('/users/whou/alsf_filbin/data/processed/logTPM_filtered.rds') ## [1:29229, 1:299] 
data <- mat
suppressMessages(library(scran))
sep <- '-'
batch <- sub(paste0(sep,'.*'),'',colnames(data))
id <- order(batch)
data <- data[,id]
batch <- batch[id]
## highly variable genes
cn <- colnames(data)
fit <- trendVar(data)
decomp <- decomposeVar(data,fit)
gs <- row.names(decomp)[decomp[,'total'] > decomp[,'tech']]
data <- data[gs,]  ###  [1:13662, 1:299]
scmd <- sapply(1:length(unique(batch)),function(sp) {
  paste0("data[gs,batch==unique(batch)[",sp,"]]")
})

# run mnn
d <- min(50,nrow(data))
ncores <- 4
cmd <- paste0("res <- fastMNN(",paste0(scmd,collapse = ","),",BPPARAM=MulticoreParam(ncores),d=",d,')')
eval(parse(text=cmd))
d <- res$corrected
row.names(d) <- cn
d <- d[,1:10]
library(umap)
set.seed(12345)
u <- umap(d)$layout
# sample effect kind of removed 
library(ggplot2)
pdf('/users/whou/alsf_filbin/data/plot/mnn_umap.pdf',width=6,height=4.5)
ggplot() + geom_point(data=data.frame(umap1=u[,1],umap2=u[,2],batch=as.factor(batch)),aes(x=umap1,y=umap2,col=batch),alpha=0.8) + theme_classic()
dev.off()

