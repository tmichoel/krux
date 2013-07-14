#A simple script to calculate FRD values 
#for each significant assocation based on results on randomized data.
#Usage:
#result.file: the result from krux based on orignial data
#random.file: the result from krux based on shuffled data
#num.randomRun: the number of runs on shuffled data
#fdr.file: the output file

result.file = '../fdr/R.output'
random.file = '../fdr/R.output.random'
num.randomRun = 10
fdr.file = 'fdr.out'

result = read.table(result.file, sep=' ', stringsAsFactors=F)
result = result[order(result$V5),]

random = read.table(random.file, sep=' ', stringsAsFactors=F)
random = sort(random[,5])
index = 1
fdr = numeric(dim(result)[1])
fdr[] =  length(random) / dim(result)[1] / num.randomRun 
for (i in seq(1,dim(result)[1])){
	#cat(i,'\n')
	while(index <= length(random) && random[index] <= result[i,5]){
		index = index + 1
	}	

        if (index > length(random)){break}

	fdr[i] = (index-1) / i / num.randomRun 	
	
}
result = cbind(result, fdr)
write.table(result, file=fdr.file, row.names=F, col.names=F, quote=F)

