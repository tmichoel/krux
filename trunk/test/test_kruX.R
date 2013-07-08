##Test script between R built-in function and matrix-based function
source('../r/kruX.R')
p.value.threshold=1
geno=read.table('kruX_testData_genotype.txt', sep='\t', stringsAsFactors=F, header=T)[,-1]
mrna=read.table('kruX_testData_expression.txt', sep='\t', stringsAsFactors=F, header=T)[,-1]
write.table(geno, file='geno.tab.tmp' ,sep='\t', quote=F)
write.table(mrna, file='mrna.tab.tmp' ,sep='\t', quote=F)

t.start = proc.time()
res1=CalculateKruskalWallisWithMatrix(geno, mrna, p.value.threshold, 1000)
cat("Total running time for R matrix-based test: ", proc.time()[3]-t.start[3],"\n")
res1=do.call(rbind, res1)
write.table(res1[order(res1$pvalue),], file='R.output', quote=F, row.names=F, col.names=F)
res1=res1[order(res1$pvalue),]

#Due to long running time, skip the comparison between 
#the results of built-in kruskal-wallis  function and 
#the result of krux. 
if (FALSE){
	res2=brutal.force.kruskal(geno, mrna, p.value.threshold)
	res2=res2[order(res2$pvalue),]

	if (T %in% (res1$pvalue-res2$pvalue)>1e-10){
  		cat('failed: R built-in function vs R matrix-based test \n')
	} else{
  		cat('passed: R built-in function vs R matrix-based test  \n')
	}
}
