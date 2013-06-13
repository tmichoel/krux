function pcut = fdrcut(preal,fdr,fdrcut)
% FDRCUT - P-value cutoff at given FDR level
% FDRCUT returns the P-value cutoff for a given FDR cutoff. Because the
% output of the "fdrvec" is not monotonic, especially at very low P-values,
% due to discrete jumps, we define the FDR cutoff as
%
%   pcut = max(preal(fdr<=fdrcut));
%
% USAGE: pcut = fdrcut(preal,fdr,fdrcut);
%
% INPUT: - preal : vector of P-values
%        - fdr : vector of FDR values, see "fdrvec" function
%        - fdrcut : FDR cutoff
%
% OUTPUT: pcut : P-value cutoff
%
% Copyright 2012-2013, Tom Michoel
%   tom.michoel@roslin.ed.ac.uk
%   http://www.roslin.ed.ac.uk/tom-michoel

pcut = max(preal(fdr<=fdrcut));