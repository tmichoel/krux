function [Iout,Jout,varargout] = cisselect(Icis,Jcis,I,J,varargin)
% CISSELECT - Select cis-pairs from kruX output
% CISSELECT filters the Kruskal-Wallis test statistics, P-values and/or
% degrees of freedom for a subset of cis-pairs from the output of kruX.
%
% USAGE: [Iout,Jout,Pcis] = cisselect(Icis,Jcis,I,J,P);
%        [Iout,Jout,Pcis,Scis] = cisselect(Icis,Jcis,I,J,P,S);
%        [Iout,Jout,Pcis,Scis,dfcis] = cisselect(Icis,Jcis,I,J,P,S,df);
%
% Copyright 2012-2013, Tom Michoel
%   tom.michoel@roslin.ed.ac.uk
%   http://www.roslin.ed.ac.uk/tom-michoel

% subset and location of cis-pairs that are in kruX output
[tf,loc] = ismember([Icis,Jcis],[I,J],'rows');
Iout = Icis(tf);
Jout = Jcis(tf);

% set output arguments
if nargout-2 == nargin-4
    varargout = cell(size(varargin));
    for k=1:length(varargin)
       varargout{k} = varargin{k}(loc(tf)); 
    end
else
    error(['Number of variable output arguments not equal to variable ' ...
         'number of input arguments']);
end