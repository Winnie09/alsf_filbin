<<<<<<< HEAD
mat = readRDS('/users/whou/alsf_filbin/data/processed/logTPM.rds') # [1:47127, 1:1152]
af <- list.files('/dcl02/hongkai/data/whou/alsf_filbin/align/')
qua <- t(sapply(af,function(f) {
  d <- readLines(paste0('/dcl02/hongkai/data/whou/alsf_filbin/align/',f,'/summary.txt'))
=======
mat = readRDS('/users/whou/alsf_filbin/data/processed/logTPM.rds')
af <- list.files('/users/whou/alsf_filbin/data/align/')
qua <- t(sapply(af,function(f) {
  d <- readLines(paste0('/users/whou/alsf_filbin/data/align/',f,'/summary.txt'))
>>>>>>> 493eb965c5992a7a722bc79e959aa16119876f2b
  c(as.numeric(sub(' .*','',d[1])),as.numeric(sub('%.*','',d[length(d)])))
}))

colnames(qua) = c('num.read','alignment_rate')
## filter cells
# alignment rate > 50%
# number of expressed genes > 500
<<<<<<< HEAD
mat = mat[,colSums(mat>0)>1000 & qua[,2]>50] 
## filter genes
mat <- mat[rowMeans(mat > 0.1, na.rm = T) > 0.01,] #  [1:27556, 1:762]
mat <- mat[!grepl('^MT-',rownames(mat)), ]
saveRDS(mat,'/users/whou/alsf_filbin/data/processed/logTPM_filtered.rds')

=======
mat = mat[,colSums(mat>0)>1500 & qua[,2]>50] ##  [1:43797, 1:299]
## filter genes
mat <- mat[rowMeans(mat > 0.1) > 0.01,] # [1:29229, 1:299]
saveRDS(mat,'/users/whou/alsf_filbin/data/processed/logTPM_filtered.rds')
>>>>>>> 493eb965c5992a7a722bc79e959aa16119876f2b
