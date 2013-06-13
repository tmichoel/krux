function f = fdrvec(preal,prand)
% FDRVEC - Convert vector of P-values to FDR values
% FDRVEC converts a vector of P-values preal to a vector of FDR values f
% based on data from one or more permutations prand. For every value p in
% preal, its FDR is defined as the ratio of the average number of random
% values below p and the number of real values below p.
%
% USAGE: f = fdrvec(preal,prand);
%
% INPUT: - preal : real pvalues
%        - prand : random pvalues, either a vector (from a single random
%          permutation) or a matrix or cell array of vectors (from a set of
%          random permutations
%
% OUTPUT: f : f(k) = FDR(preal(k))
%
% Copyright 2012-2013, Tom Michoel
%   tom.michoel@roslin.ed.ac.uk
%   http://www.roslin.ed.ac.uk/tom-michoel

N = length(preal);

if ~iscell(prand) && size(prand,2)==1
    % work on sorted data
    [p1,I] = sort(preal,'ascend');
    p2 = sort(prand,'ascend');
    % sort jointly
    [~,J] = sort([p1;p2],'ascend');
    % wherever J>N we have a random one inserted in the real order
    nrand = cumsum(J>N);
    % back to index for real values
    nreal = (1:N)'; % number of real values below each entry of p1
    nrand = nrand(J<=N); % number of random values below each entry of p1
    % FDR value
    f = nrand./nreal;
    % permute back to original index
    [~,Iinv] = sort(I,'ascend'); % inverse permutation
    f = f(Iinv);
elseif iscell(prand) % random data in cell array
    fset = zeros(N,length(prand));
    for k=1:length(prand)
       fset(:,k) = fdrvec(preal,prand{k}); 
    end
    f = mean(fset,2);
else % random data as column vectors in a matrix
    fset = zeros(N,size(prand,2));
    for k=1:size(prand,2)
       fset(:,k) = fdrvec(preal,prand(:,k)); 
    end
    f = mean(fset,2);
end