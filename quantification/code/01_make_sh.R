allfd = list.files('/dcl02/hongkai/data/whou/filbin-pediatric-high-grade-glioma/fastq')
for (fd in allfd){
  print(fd)
  allfn = list.files(paste0('/dcl02/hongkai/data/whou/filbin-pediatric-high-grade-glioma/fastq/',fd,'/'))
  existfn = list.files(paste0('/users/whou/alsf_filbin/data/sample/',fd))
  allfn = setdiff(allfn, existfn)
  allfn = gsub('_R1.fastq.gz','',allfn)
  allfn = unique(gsub('_R2.fastq.gz','',allfn))
  print(length(allfn))
  for (fn in allfn){
    file = paste0(fn,'.sh')
    sink(paste0('/users/whou/alsf_filbin/shfiles/',file))
    writeLines(paste0('mkdir /users/whou/alsf_flibin/data/align/',fn))
    writeLines(paste0('cd /users/whou/alsf_flibin/data/align/',fn))
    writeLines(paste0('/users/whou/software/hisat2/hisat2-2.1.0/hisat2 -x /users/whou/software/hisat2/genome/hg19/genome -1 /dcl02/hongkai/data/whou/filbin-pediatric-high-grade-glioma/fastq/', fd, '/', fn, '_R1.fastq.gz -2 /dcl02/hongkai/data/whou/filbin-pediatric-high-grade-glioma/fastq/', fd, '/', fn, '_R2.fastq.gz -S data.sam 2> summary.txt'))
    writeLines('/jhpce/shared/community/core/samtools/1.1/bin/samtools view -bS data.sam > data_unsorted.bam')
    writeLines('/jhpce/shared/community/core/samtools/1.1/bin/samtools sort data_unsorted.bam data')
    writeLines('rm data.sam')
    writeLines('rm data_unsorted.bam')
    writeLines('/users/whou/software/stringtie/stringtie-2.0.3.Linux_x86_64/stringtie data.bam -G /users/whou/software/stringtie/gencode.v19.annotation.gtf -e -A genes.out -o data.gtf')
    sink()
  }
}

