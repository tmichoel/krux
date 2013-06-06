"""Matrix-based implementation for Kruskal-wallis test
Tested on Python3.2 with numpy1.7 and scipy0.12
"""

import numpy 
import scipy.stats
#import IPython
import pdb
import subprocess
import os 
import time

def calculateKruskalWallisWithMatrix(geno, 
                                     mrna,
                                     pValueThre = 1e-3,
                                     numTranscript = 1000,
                                     shuffledMrna = False):
  """Compute Kruskal-Wallis test p-values for pairs of transcripts and SNPs with
  matrix operations.
  Args:
             geno: Matrix for SNPs where each column are SNPs of a particular
                   sample and each row are the values of a particular SNP 
                   cross all samples. 
             mrna: Matrix for expressiion values of transcrtips. The format 
                   is defined same as genotype matrix.
       pValueThre: Only pairs of SNP-transcript  with p-values less than the 
                   threshold are included in the calculation
    numTranscript: The number of transcripts to be calculated in one round.
                   The parameter is determined by available  memory size  
     shuffledMrna: If T, mrna is shuffled, so the result is based on 
                   perturbation data. Default value is F
  Returns: 
     Ordered transcript-SNP pairs by their p-values
  
  version 1.0 May 2, 2013 Initial code
  """

  def preCalculate(snpsID, mrnaStart, mrnaEnd):
    """This function calculates freedome for chi-square test, numbers of SNPs 
    for no-missing  transcripts, and threshold for the corresponding p-value
    Note that the number of SNPs with a particular value
    should be corrected for missing values in transcripts
    """
    
    #mark mrna missing value 
    mrnaMissing = numpy.zeros(mrna[mrnaStart:mrnaEnd].shape)
    mrnaMissing[numpy.where(numpy.isnan(mrna[mrnaStart:mrnaEnd]))] = 1

    #Calculate the number of missing transcripts for each SNP with a value from 0,1,2 
    snpsCorrect = genoRank[:,snpsID[0],:].dot(mrnaMissing.transpose())

    #The degree of freedom of chi-square test between a pair of SNP-transcript is the number 
    #of different values of the SNP in all samples with non-missing value of the transcript
    #N is the effective number of samples for each test
    #chiSquareFD is the freedom for each test
    #Broadcasting snps[snpsID] to match snpsCorrect 
    N = snps[snpsID].transpose()[:,:,numpy.newaxis]-snpsCorrect 
    chiSquareFD = numpy.apply_along_axis(numpy.sum, 0, N>0) - 1 
    N = numpy.apply_along_axis(numpy.sum, 0, N)

    #Threshold for the chi-square test between each pair given the p-value
    chi2Thre = (scipy.stats.chi2.ppf(1-pValueThre, chiSquareFD) +  
               3 * (N + 1)) * (N * (N + 1)) / 12 
    chi2Thre[numpy.where(chiSquareFD<1)] = numpy.iinfo(numpy.int32).max

    return dict(mrna=(mrnaStart, mrnaEnd), snpsCorrect=snpsCorrect,
                N=N, chiSquareFD=chiSquareFD, chi2Thre=chi2Thre)

  def calculateNoMissing(snpsID, preResult):
    """This function handles SNPs without missing value.
    """

    #Calculate ranks for selected transcripts. NA's rank is replaced by 0
    #Note that masked array is also an option
    mrnaSelected = mrna[preResult['mrna'][0]:preResult['mrna'][1]]
    expressionRank = numpy.apply_along_axis(scipy.stats.rankdata,1,
                      mrnaSelected)
    expressionRank[numpy.where(numpy.isnan(mrnaSelected))] = 0
    
    #calulate rank statistics for different genotype values(0,1,2) and then sum
    testStat = numpy.power(numpy.dot(genoRank[:,snpsID[0],:], expressionRank.transpose()),2) / \
               (snps[snpsID].transpose()[:,:,numpy.newaxis] - preResult["snpsCorrect"])
    testStat[numpy.where(numpy.isnan(testStat))] = 0
    testStat = numpy.apply_along_axis(numpy.sum, 0,testStat) 
    preResult['testStat'] = testStat

    return

  def calculateMissing(snpsID, preResult):
    """This function handles SNPs with missing value.
    Note that only one transcript is processed in each round 
    when working on SNPs with missing values.
    """
    
    #Order SNPs according to the transcript order (nan is last) 
    #advanced indexing of genoRank
    mrnaSelected = mrna[preResult['mrna'][0]:preResult['mrna'][1]]
    mrnaOrder = numpy.argsort(numpy.ravel(mrnaSelected)) 
    genoRankOrder = genoRank[:,numpy.ravel(snpsID),:][:,:,mrnaOrder]

    #Calculate correction due to missing value in SNPs 
    sumCor = geno[numpy.ravel(snpsID),:][:,mrnaOrder]
    cumSumCor = numpy.where(numpy.isnan(sumCor), 1, 0)
    cumSumCor = numpy.cumsum(cumSumCor, axis=1)
 
    #Since missing value in the transcript are represented by 
    #zeor in order, the corresponding part in culcumSumCor is replaced by zeor, too. 
    cumSumCor = cumSumCor * (-numpy.sort(numpy.where(
                numpy.isnan(mrnaSelected), 0, -1)))
    orderStat = numpy.arange(1,1+mrnaOrder.size) * \
                (-numpy.sort(numpy.where(numpy.isnan(mrnaSelected), 0, -1)))

    #calulate rank statistics for different genotype values(0,1,2) and then sum 
    testStat= numpy.power(
                 (numpy.dot(genoRankOrder, orderStat.transpose())-
                 numpy.sum(genoRankOrder * cumSumCor[numpy.newaxis,:,:], axis=2)
                 [:,:,numpy.newaxis]), 2) / \
              (snps[snpsID].transpose()[:,:,numpy.newaxis] - 
               preResult["snpsCorrect"])
    testStat[numpy.where(numpy.isnan(testStat))] = 0
    testStat = numpy.apply_along_axis(numpy.sum, 0,testStat)
    preResult['testStat'] = testStat

    return

  def postCalculate(snpsID, preResult):
    """This function caculates P-vauels for pairs with statistics more than 
    the threshold and save the result to output list
    """

    #find the index of pairs with pvalue less then threshold
    #build record array to save information of selected pairs
    #Jun 2, 2013 Bug fixed >= instead of >
    pairSelected = numpy.where(preResult['testStat'] >= preResult['chi2Thre'])
    SNPs = numpy.ravel(snpsID)[pairSelected[0]]
    genes = numpy.arange(preResult['mrna'][0], preResult['mrna'][1])[pairSelected[1]]
    fd = preResult['chiSquareFD'][pairSelected]
    chi = (preResult['testStat'][pairSelected] * 12 / 
          (preResult['N'][pairSelected] * (preResult['N'][pairSelected] + 1))
          - 3 * (preResult['N'][pairSelected] + 1))
    pValue = scipy.stats.chi2.sf(chi, fd) #use sf to get more accrate value

    #return a structured array
    result = numpy.core.records.fromarrays([SNPs, genes, fd, chi, pValue],
             names='SNP,gene,fd,chi,pvalue')
    return result

  #list for the result
  output = list()

  #seperate SNP indexs by if having missing values
  snpsWithMissing = numpy.where(numpy.any(numpy.isnan(geno), 1)) 
  snpsNoMissing = numpy.where(numpy.all(numpy.logical_not(numpy.isnan(geno)), 1))
  
  #Genotype position index matrix  with a 3 dimensional array
  genoRank = numpy.array([numpy.where(geno==i,1,0)  
             for i in numpy.unique(numpy.ma.masked_invalid(geno).compressed())]) 

  #The number of SNPs with a particular value
  snps = numpy.sum(genoRank, axis=2).transpose()
  
  #Calculation for SNPs without missing values
  if snpsNoMissing[0].size>0:
    mrnaStart = list(range(0,mrna.shape[0],numTranscript))
    mrnaEnd   = mrnaStart[1:] + [mrna.shape[0]] 
    for ids in zip(mrnaStart, mrnaEnd):
      #Process a slice of  mrna
      preResult = preCalculate(snpsNoMissing, ids[0], ids[1])
      calculateNoMissing(snpsNoMissing, preResult) 
      result = postCalculate(snpsNoMissing, preResult) 
      output.append(result) if result.size>0 else None
  
  #Calculation for SNPs with missing values
  if snpsWithMissing[0].size>0:
    for mrnaIX in range(0, mrna.shape[0]):
      #Process a single mrna
      preResult = preCalculate(snpsWithMissing, mrnaIX, mrnaIX+1)
      calculateMissing(snpsWithMissing, preResult)
      result = postCalculate(snpsWithMissing, preResult)
      output.append(result) if result.size>0 else None

  if (output):#Bug
    output = numpy.concatenate(output)
    output = output[numpy.argsort(output['pvalue'])]
  
  return(output)


if __name__ == '__main__':
  for i in range(10):
    print(str(i))
