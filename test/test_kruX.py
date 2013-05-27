import numpy
import subprocess
import time
import os
import sys
sys.path.append('..')
import kruX

def test(numSample, numSNP, numMRNA, perMissing, pThreshold):
  """ Test the results from R's build-in kruskal wallis function, R's matrix-based
  kw test, and Python's matrix-based kw test
  directory.  Can be done in a nice way with rpy
  """
  subprocess.call(['Rscript','test_kruX.R',str(numSample),
                   str(numSNP), str(numMRNA), str(perMissing), str(pThreshold)])
  outputFromR = numpy.genfromtxt('R.output',
                dtype=[('SNP','i4'),('gene','i4'),('fd','i4'),
                       ('chi','f4'),('pvalue','f4') ] )
  geno = numpy.genfromtxt('geno.tab.tmp', skip_header=1, missing_values='NA')
  geno = geno[:,1:]
  mrna = numpy.genfromtxt('mrna.tab.tmp', skip_header=1, missing_values='NA')
  mrna = mrna[:,1:] 
  start = time.clock()
  outputFromPython = kruX.calculateKruskalWallisWithMatrix(geno, mrna, 
                     pValueThre = pThreshold)
  endTime = time.clock() - start
  print('Total running time for Python matrix-based test: %.2f'%(endTime))
  numpy.savetxt('python.output', outputFromPython, fmt='%i\t%i\t%i\t%.7f\t%.7f',delimiter='\t')
  if (numpy.all(outputFromR['pvalue']-outputFromPython['pvalue']<1e-7)):
    print('Passed: equal between R matrix-based  and Python matrix-based test')
  else:  
    print('Failed: not equal between R matrix-based  and Python matrix-based test')

def test_python_performance(numSample, numSNP, numMRNA, perMissing, pThreshold):
  """ Test the results from R's build-in kruskal wallis function, R's matrix-based
  kw test, and Python's matrix-based kw test
  directory.  Can be done in a nice way with rpy
  """
  geno = numpy.genfromtxt('geno.tab.tmp', skip_header=1, missing_values='NA')
  geno = geno[:,1:]
  mrna = numpy.genfromtxt('mrna.tab.tmp', skip_header=1, missing_values='NA')
  mrna = mrna[:,1:]
  start = time.clock()
  outputFromPython = kruX.calculateKruskalWallisWithMatrix(geno, mrna,
                     pValueThre = pThreshold)
  endTime = time.clock() - start
  print('Total running time for Python matrix-based test: %.2f'%(endTime))
  numpy.savetxt('python.output', outputFromPython, fmt='%i\t%i\t%i\t%.7f\t%.7f',delimiter='\t')

if __name__ == '__main__':
  for i in range(1):
    print(str(i))
    test(100, 100, 100, 0.2, 0.8)
