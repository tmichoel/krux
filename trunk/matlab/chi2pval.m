function p = chi2pval(x,v)
% CHI2PVAL - Compute chi2 P-value
% CHI2PVAL computes the chi2 P-value using the expression found <a
% href="http://www.mathworks.co.uk/support/solutions/en/data/1-DPEPZF/index.html">here</a>,
% instead of using 1-ch2cdf(x,v), to get higher accuracy at very low
% P-values. This function also works without the Statistics Toolbox.

p = gammainc(x/2,v/2,'upper');

