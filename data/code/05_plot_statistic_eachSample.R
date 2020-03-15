mat = readRDS('/users/whou/alsf_filbin/data/processed/logTPM.rds')
sep <- '-'
batch <- sub(paste0(sep,'.*'),'',colnames(mat))

af <- list.files('/users/whou/alsf_filbin/data/align/')
qua <- t(sapply(af,function(f) {
  d <- readLines(paste0('/users/whou/alsf_filbin/data/align/',f,'/summary.txt'))
  c(as.numeric(sub(' .*','',d[1])),as.numeric(sub('%.*','',d[length(d)])))
}))
colnames(qua) = c('num.read','alignment_rate')

pdf('/users/whou/alsf_filbin/data/plot/statistic_eachSample.pdf',width=10,height=15)
par(mfrow=c(6,4))
for (i in unique(batch)){
  hist(colSums(mat)[batch==i], col='grey', breaks=50, main=i,xlab='sum.log2TPM',ylab='num.cells')
  hist(colSums(mat>0)[batch==i],col='grey', breaks=50, main=i,xlab='num.genes.detected',ylab='num.cells')
  hist(qua[batch==i,'num.read'],col='grey',breaks=50,main=i,xlab='num.reads',ylab='num.cells')
  hist(qua[batch==i,'alignment_rate'],col='grey',breaks=50,main=i,xlab='alignment.rate',ylab='num.cells')
}
dev.off()

