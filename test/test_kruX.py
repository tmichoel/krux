import numpy
import subprocess
import time
import os
import sys
sys.path.append('../python/')
import kruX

def test():
  """ Test the results from R's build-in kruskal wallis function, R's matrix-based
  kw test, and Python's matrix-based kw test
  directory.  Can be done in a nice way with rpy
  """
  subprocess.call(['Rscript','test_kruX.R'])
  outputFromR = numpy.genfromtxt('R.output',
                dtype=[('SNP','i4'),('gene','i4'),('fd','i4'),
                       ('chi','f4'),('pvalue','f4') ] )
  geno = numpy.genfromtxt('geno.tab.tmp', skip_header=1, missing_values='NA')
  geno = geno[:,1:]
  mrna = numpy.genfromtxt('mrna.tab.tmp', skip_header=1, missing_values='NA')
  mrna = mrna[:,1:] 
  start = time.clock()
  outputFromPython = kruX.calculateKruskalWallisWithMatrix(geno, mrna, 
                     pValueThre = 1, numTranscript=1000)
  endTime = time.clock() - start
  print('Total running time for Python matrix-based test: %.2f'%(endTime))
  numpy.savetxt('python.output', outputFromPython, fmt='%i\t%i\t%i\t%.17e\t%.17e',delimiter='\t')
  if (numpy.all(outputFromR['pvalue']-outputFromPython['pvalue']<1e-7)):
    print('Passed: equal between R matrix-based  and Python matrix-based test')
  else:  
    print('Failed: not equal between R matrix-based  and Python matrix-based test')

def test_python_performance():
  """Python's matrix-based kw test running time
  """
  geno = numpy.genfromtxt('geno.tab.tmp', skip_header=1, missing_values='NA')
  geno = geno[:,1:]
  mrna = numpy.genfromtxt('mrna.tab.tmp', skip_header=1, missing_values='NA')
  mrna = mrna[:,1:]
  start = time.clock()
  outputFromPython = kruX.calculateKruskalWallisWithMatrix(geno, mrna,
                     pValueThre = 1)
  endTime = time.clock() - start
  print('Total running time for Python matrix-based test: %.2f'%(endTime))
  numpy.savetxt('python.output', outputFromPython, fmt='%i\t%i\t%i\t%.7f\t%.7f',delimiter='\t')

if __name__ == '__main__':
  test()
  #test_python_performance()
