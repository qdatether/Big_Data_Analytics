function  xlagM = multilagmatrix(xM,responseindex,ordersV,indexV)
% multilagmatrix builds the set of explanatory variables for the dynamic
% regression model.
% INPUTS
% - xM          : the matrix of $K$ time series (variables in columns)
% - responseindex : the index of the response variable in {1,...,K}
% - ordersV     : vector of size 1xK of the maximum order for each of the K
%                 variables.
% - indexV      : the vector of size 1 x K*pmax of zeros and ones
%                 e.g. if the component in position 2*pmax+3 is one, the
%                 third variable, lag 3, X3(t-3), is selected.
% OUTPUTS
% - xlagM       : the matrix of all explanatory lagged variables in the
%                 DR model. The sequence of the lagged variables in 'lagM'
%                 is determined by indexV.
[n,K]=size(xM);
pmax=size(indexV,2)/K;
xtempM=NaN(n,K*pmax);
for iK=1:K
    xtempM(:,((iK-1)*pmax+1):(ordersV(iK)+((iK-1)*pmax)))=...
        lagmatrix(xM(:,iK),[1:ordersV(iK)]);
end
xlagM = [xM(:,responseindex) xtempM(:,find(indexV==1))];
xlagM = xlagM(max(ordersV)+1:end,:);
