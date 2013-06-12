function varargout = kruX(varargin)
% kruX - Matrix algorithm to compute Kruskal-Wallis test statistics
% kruX computes the Kruskal-Wallis test statistic for every pair of
% rows in D (matrix of expression values; missing data indicated with NaN)
% and C (matrix of genotypes taking values 0,1,2,...; missing data can have
% value -1 or NaN).
%
% USAGE: [I,J,P,S,df] = kruX(D,C,Pcut);
%
%       Retains all pairs (I,J) where P<=Pcut; S is the corresponding 
%       Kruskal-Wallis test statistic and df the degrees of freedom for the
%       chi2-test (number of genotype values - 1 for each row).
%
%       [I,J,S,df] = kruX(D,C,Scut,flag);
%
%       If you don't have Statistics Toolbox, give dummy input variable
%       flag = 1 and provide vector of cutoffs for test statistic (one for
%       each possible degree of freedom). All pairs (I,J) where S exceeds
%       the appropriate cutoff are returned, without P-values.
%
% DEPENDENCY: kruX with 3 input arguments uses the Statistics Toolbox to
%             compute tied ranks and convert test statistics to P-values;
%             with 4 input arguments there are no external dependencies,
%             but no correction for ties is made and no P-values are
%             reported.
%
% Copyright 2012-2013, Tom Michoel
%   tom.michoel@roslin.ed.ac.uk
%   http://www.roslin.ed.ac.uk/tom-michoel

switch nargin
    case 3
        D = varargin{1};
        C = varargin{2};
        Pcut = varargin{3};
        flag = 0;
    case 4
        D = varargin{1};
        C = varargin{2};
        Scut = varargin{3};
        flag = 1;
    otherwise
        error('Wrong number of input arguments');
end

if flag==0 % we have statistics toolbox
    [R, tieadj] = tiedrank(D'); % replace values in D by their rank in each row
    R = R';
    R(isnan(D)) = 0;
else % without statistics toolbox we don't take tied ranks into account
    R = rankrows(D);
end

K = size(D,2); % number of samples

missD = nnz(R==0); % number of missing expression values
if missD>0
    Z = double(R~=0); % indicator for non-missing values
end

% number of genotype values
nC = double(max(C(:)));
% cutoffs for test statistic for different degrees of freedom
if flag==0
    Scut = chi2inv(1-Pcut,1:nC);
end

% create genotype group index vectors
Ic = cell(nC+1,1);
for k=0:nC
   Ic{k+1} = sparse(C==k); 
end
% same for missing values
In = C==-1 | isnan(C);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Genotypes without missing values %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ixc = find(sum(In,2)==0);

% there may be no such markers ...
if ~isempty(ixc)

    % number of elements in each group for all pairs of rows
    Ne = cell(nC+1,1);
    if missD==0
        for k=0:nC
            Ne{k+1} = repmat(sum(Ic{k+1}(ixc,:),2)', [size(R,1),1]);
        end
    else
        for k=0:nC
            Ne{k+1} = Z*Ic{k+1}(ixc,:)';
        end
    end
    
    % sum of (mean rank)^2 x number of elements in each group for all pairs of
    % rows
    S = zeros(size(D,1),length(ixc));
    for k=0:nC
        mR = ((R*Ic{k+1}(ixc,:)').^2)./Ne{k+1};
        % any NaN result from 0/0 for empty groups, set them to 0 to remove from sum
        mR(isnan(mR)) = 0;
        % sum
        S = S + mR;
    end
    
    % Kruskal-Wallis test statistic
    if missD==0
        S = 12*S./(K.*(K+1)) - 3*(K+1);
        % With Statistics Toolbox we can adjust for ties
        if flag==0
            adjust = 1 - 2*tieadj'./(K^3-K);
            S = S./repmat(adjust,[1 size(S,2)]);
        end
    else
        % number of xpr samples for each pair
        Kmat = repmat(sum(Z,2),[1 length(ixc)]);
        S = 12*S./(Kmat.*(Kmat+1)) - 3*(Kmat+1);
        % With Statistics Toolbox we can adjust for ties
        if flag==0
            adjust = 1 - 2*repmat(tieadj',[1 size(S,2)])./(Kmat.^3-Kmat);
            S = S./adjust;
        end
    end
    
    % degrees of freedom is number of groups minus 1
    df = -ones(size(S,1),size(S,2));
    for k=0:nC
        df = df + spones(Ne{k+1});
    end
    
    % find pairs above cutoff
    bool = S>=Scut(1) & df==1;
    for k=2:nC
        bool = bool | (S>=Scut(k) & df==k);
    end
    [I,J] = find(bool);
    ind = sub2ind(size(S),I,J);
    S = S(ind);
    df = df(ind);
    J = ixc(J); % back to original index
    if size(I,1)==1
        I = I';
        S = S';
        df = df';
    end

else
    % make empty arrays
    I = [];
    J = [];
    S = [];
    df = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Genotypes with missing values %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ixn = find(sum(In,2)>0);

% unique missing data patterns
[Un, ~, Bn] = unique(In(ixn,:),'rows');

% do calculation for each unique row
for k=1:size(Un,1)
    % new markers
    mloc = ixn(Bn==k);
    Cnew =  C(mloc, Un(k,:)==0);
    % new data
    Dnew = D(:, Un(k,:)==0);
    % recursive call
    if nnz(isnan(Cnew) | Cnew==-1)>0
        error('this cannot be true')
    end
    if flag==0
        [Inew,Jnew,~,Snew,dfnew] = kruX(Dnew,Cnew,Pcut);
    else
        [Inew,Jnew,Snew,dfnew] = kruX(Dnew,Cnew,Scut,1);
    end
    % change back to original index
    Jnew = mloc(Jnew);
    % append results
    I = [I;Inew];
    J = [J;Jnew];
    S = [S;Snew];
    df = [df;dfnew];
end

%%%%%%%%%%%%%%
%%% Output %%%
%%%%%%%%%%%%%%

% sort nicely
[I,t] = sort(I);
J = J(t);
S = S(t);
df = df(t);

switch nargout
    case 4
        varargout{1} = I;
        varargout{2} = J;
        varargout{3} = S;
        varargout{4} = df;
    case 5
        varargout{1} = I;
        varargout{2} = J;
        if flag==0
            % compute P-values if we have statistics toolbox
            varargout{3} = 1-chi2cdf(S,df);
        else
            warning('Cannot output P-values without Statisics Toolbox')
            varargout{3} = nan(size(S));
        end
        varargout{4} = S;
        varargout{5} = df;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Function to rank data if Statistics Toolbox is unavailable %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [R,T] = rankrows(M)
% RANKROWS - R(i,j) is rank of M(i,j) in row i
% RANKROWS converts a matrix M to a matrix of ranks R such that R(i,j) is
% the rank of M(i,j) in row i. 
%
% USAGE: [R,T] = rankrows(M);
%        R : matrix of row ranks
%        T : matrix of permutations to sort rows
%
% WARNING: rankrows makes no correction for ties
%
% Copyright 2012, Tom Michoel
%   tom.michoel@roslin.ed.ac.uk
%   http://www.roslin.ed.ac.uk/tom-michoel

[~,T] = sort(M,2);
[~,R] = sort(T,2);
% missing data (NaN's) are sorted last and get highest rank;
% setting them to zero gives correct ranks for future calculations
R(isnan(M)) = 0;


