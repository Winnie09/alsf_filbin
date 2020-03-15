library(data.table)
<<<<<<< HEAD
af <- list.files('/dcl02/hongkai/data/whou/alsf_filbin/align')
res <- sapply(af,function(f) {
  print(f)
  d <- fread(paste0('/dcl02/hongkai/data/whou/alsf_filbin/align/',f,'/genes.out'),data.table = F)
=======
af <- list.files('/users/whou/alsf_filbin/data/align')
res <- sapply(af,function(f) {
  print(f)
  d <- fread(paste0('/users/whou/alsf_filbin/data/align/',f,'/genes.out'),data.table = F)
>>>>>>> 493eb965c5992a7a722bc79e959aa16119876f2b
  tmp <- d[,'TPM']
  names(tmp) <- paste0(d[,2],'_',d[,1])
  tmp[tmp > 0] ###
})

cn <- names(res)
gn <- unique(unlist(sapply(res,names)))
mat <- matrix(0,nrow=length(gn),ncol=length(cn),dimnames = list(gn,cn))
for (i in names(res)) {
  mat[names(res[[i]]),i] <- res[[i]]
}
saveRDS(mat,file='/users/whou/alsf_filbin/data/processed/TPM.rds')
mat <- log2(mat + 1)
saveRDS(mat,file='/users/whou/alsf_filbin/data/processed/logTPM.rds')
<<<<<<< HEAD
=======

>>>>>>> 493eb965c5992a7a722bc79e959aa16119876f2b
