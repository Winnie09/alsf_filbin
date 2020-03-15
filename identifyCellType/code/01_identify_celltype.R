source('/users/whou/alsf_filbin/R/processCountMatrix.R')
af <- c('MouseAdultDAMoleculeCounts.rds','iPSC_fulldataset.rds','Human_Embryo_fulldataset.rds','ES_fulldataset.rds')##'Mouse_Embryo_fulldataset.rds'
d <- lapply(af, function(f){
  print(f)
  m = readRDS(paste0('/users/whou/alsf_filbin/data/atlases/LaManno_dataset/proc/expr/',f))
  a = readRDS(paste0('/users/whou/alsf_filbin/data/atlases/LaManno_dataset/proc/meta/',f))
  tmp <- sapply(unique(a$celltype), function(i) rowSums(m[, a$celltype==i,drop=F]))   
  processCountMatrix(tmp, platform='fluidigm', log2NormCount = TRUE, useHomologGene = TRUE, species=ifelse(grepl('Human',f),'human','mouse'), removeDupGeneName = TRUE) 
})
names(d) = af

gn <- sapply(d,row.names)
gn <- table(unlist(gn))
gn <- names(gn)[gn==length(d)]
for (i in 1:length(d)) d[[i]] <- d[[i]][gn,,drop=F]
data = do.call(cbind, d)


m = read.table('/users/whou/alsf_filbin/data/Kriegstein_Human_cortex/exprMatrix.tsv',as.is=T)
m = as.matrix(m)
rownames(m) = m[,1]
m = m[,-1]
colnames(m) = m[1,]
m = m[-1,]
dn = dimnames(m)
m = matrix(as.numeric(m),nrow=nrow(m))
dimnames(m) = dn

a = read.table('/users/whou/alsf_filbin/data/atlases/Kriegstein_Human_cortex/meta.tsv',as.is=T, fill=T, header=T)
a = a[complete.cases(a),]
ct = a[match(colnames(m), a$Cell),2]
names(ct) =colnames(m)
ct = ct[!grepl('Sample',ct)]
m = m[,names(ct)]
tmp = sapply(unique(ct), function(i) rowSums(m[, ct == i, drop=F]))
data2 <- processCountMatrix(tmp, platform='fluidigm', log2NormCount = TRUE, useHomologGene = TRUE, species='human', removeDupGeneName = TRUE) 

gn = intersect(rownames(data), rownames(data2))
ref.data = cbind(data[gn,], data2[gn,])
dn = dimnames(ref.data)
library(preprocessCore)
ref.data = normalize.quantiles(ref.data)
dimnames(ref.data) = dn

test = readRDS('/users/whou/alsf_filbin/data/processed/logTPM_filtered.rds')
rownames(test) = sub('_.*','',rownames(test))
gn = intersect(rownames(test), rownames(ref.data))
test = test[gn,]
ref.data = ref.data[gn,]
ref.data.bak = ref.data
test.bak = test
# test = test[,1:100]

library(SingleR)
pred = SingleR(test = test, ref = ref.data, labels = colnames(ref.data))
pred2 = do.call(cbind,as.matrix(pred)@listData)
pred2 = pred@listData$scores 
rownames(pred2) = colnames(test)
donor = sub('-.*','',colnames(test))

labels = names(sort(table(pred$labels),decreasing = T)[1:20])
pdf('/users/whou/alsf_filbin/identifyCellType/plot/hm.pdf',width=12,height=12.5)
plotScoreHeatmap(pred, clusters=donor)
dev.off()
