countfunc <- function(mat, platform='fluidigm', normalize = TRUE,species='human', log.transform.only = FALSE) {
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

