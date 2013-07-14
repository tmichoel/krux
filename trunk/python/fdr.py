"""A simple script to calculate FRD values 
for each significant assocation based on results on randomized data.
Usage:
result_file: the result from krux based on orignial data
random_file: the result from krux based on shuffled data
num_randomRun: the number of runs on shuffled data
fdr_file: the output file
"""

import numpy
import numpy.lib.recfunctions

def cal_fdr(p_value, p_value_shuffled, num_run):

  fdr = numpy.empty(p_value.shape)
  fdr[:] = p_value_shuffled.size / p_value.size / num_run   
  index = 0
  for i in range(p_value.shape[0]):
    while (index < p_value_shuffled.shape[0]) and \
          (p_value_shuffled[index] <= p_value[i]):
      index += 1  
    if index >= p_value_shuffled.shape[0]:
      break
    fdr[i] = index / (i+1) / num_run

  return fdr

if __name__ == '__main__':
  result_file = '../fdr/R.output'
  random_file = '../fdr/R.output.random'  
  num_randomRun = 10
  fdr_fie = 'fdr.out'

  result = numpy.genfromtxt(result_file,
           dtype=[('SNP','S10'),('gene','S10'),('fd','i4'),
           ('chi','f4'),('pvalue','f4') ] )
  result = result[numpy.argsort(result['pvalue'])] 
  random_result = numpy.genfromtxt(random_file,
                   dtype=[('SNP','S10'),('gene','S10'),('fd','i4'),
                   ('chi','f4'),('pvalue','f4') ] )
  random_result = random_result[numpy.argsort(random_result['pvalue'])] 
  fdr_value = cal_fdr(result['pvalue'], random_result['pvalue'], num_randomRun)
  result= numpy.lib.recfunctions.append_fields(result, names='fdr',data=fdr_value,dtypes='f')
  numpy.savetxt(fdr_fie,result,
                fmt='%s\t%s\t%i\t%.5e\t%.5e\t%.5e',delimiter='\t')
