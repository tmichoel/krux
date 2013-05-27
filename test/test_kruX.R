##Test script between R built-in function and matrix-based function
source('../kruX.R')
args <- commandArgs(trailingOnly = TRUE)
num.sample=as.integer(args[1])
num.snp=as.integer(args[2])
num.mrna=as.integer(args[3])
per.missing = as.numeric(args[4])
p.value.threshold = as.numeric(args[5])

cat(args[1], ' samples \n', 
    args[2],     ' SNPs    \n',
    args[3],    ' transcripts \n' ,
    args[4] , ' missing values \n' ,
    args[5] , ' pvalue threshold \n')
geno=matrix(sample(c(0,1,2, NA), size=num.sample*num.snp, 
                  replace=T, prob=c(0.32,0.32,0.32,0.04)),
                 nrow=num.snp, ncol=num.sample,
                 dimnames=list(seq(0,num.snp-1), seq(0,num.sample-1))) 

mrna=matrix(rnorm(num.mrna*num.sample),
                  nrow=num.mrna, ncol=num.sample,
                  dimnames=list(seq(0,num.mrna-1), seq(0,num.sample-1)))
mrna[cbind(sample(dim(mrna)[1],round(num.mrna*num.sample*per.missing), replace=T), 
     sample(dim(mrna)[2],round(num.mrna*num.sample*per.missing), replace=T))] <- NA

write.table(geno, file='geno.tab.tmp' ,sep='\t', quote=F)
write.table(mrna, file='mrna.tab.tmp' ,sep='\t', quote=F)

t.start = proc.time()
res1=CalculateKruskalWallisWithMatrix(geno, mrna, p.value.threshold)
res1=do.call(rbind, res1)
cat("Total running time for R matrix-based test: ", proc.time()[3]-t.start[3],"\n")
write.table(res1[order(res1$pvalue),], file='R.output', quote=F, row.names=F, col.names=F)
res1=res1[order(res1$pvalue),]

res2=brutal.force.kruskal(geno, mrna, p.value.threshold)
res2=res2[order(res2$pvalue),]

if (T %in% (res1$pvalue-res2$pvalue)>1e-10){
  cat('failed: R built-in function vs R matrix-based test \n')
} else{
  cat('passed: R built-in function vs R matrix-based test  \n')
}
