function [indexV,maxorder,MSEval]=mBTS(xM,responseindex,pmax)
% mBTS implements the algorithm of modified Backward-in-Time Selection
% (mBTS), which builds a dynamic regression (DR) model selecting for
% inclusion in the model only the most relevant lagged variables to the
% response. The original BTS algorithm is presented in
% Vlachos I and Kugiumtzis D (2013), "Backward-in-time selection of the
% order of dynamic regression prediction model", Journal of Forecasting,
% 32: 685–701.
% The difference of mBTS to BTS is that BTS selects all the lags for one
% variable up to the largest lag found significantly explanatory, while
% mBTS selects only the lags found significantly explanatory (skipping
% smaller lags of the same variable if not found explanatory).
% The algorithm takes as input a matrix 'xM' of K time series (variables
% in columns), the index of the response variable ('responseindex' takes
% a value from 1 to K), and a maximum order 'pmax' to search for lags.
% INPUTS
% - xM          : the matrix of time series (the K columns are for the
%                 variables and the rows for the observations)
% - responseindex : the index of the response variable in {1,...,K}
% - pmax        : the maximum order to search for lags
% OUTPUTS
% - indexV      : the vector of size 1 x K*pmax of zeros and ones
%                 e.g. if the component in position 2*K+3 is one, the
%                 second variable, lag 3, X2(t-3), is selected.
% - maxorder    : gives the maximum order of the dynamic regression model,
%                 i.e. the maximum lag selected.
% - MSEval      : the mean of squared errors of the fitted values (the
%                 variance of the residuals).
%
%     Copyright (C) 2021 Dimitris Kugiumtzis
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%=========================================================================
% Reference : E. Siggiridou, D. Kugiumtzis, "Granger Causality in 
% Multi-variate Time Series using a Time Ordered Restricted Vector 
% Autoregressive Model", IEEE Transactions on Signal Processing, Vol 64(7), 
% pp 1759-1773, 2016
% Link      : http://users.auth.gr/dkugiu/
%=========================================================================

[n,K]=size(xM);
for d=1:K                       % remove the mean of each time series
    xM(:,d)=xM(:,d)-mean(xM(:,d));
end
indexV=zeros(1,K*pmax); % Initially, no lagged variable is selected, the 
        % first pmax positions for variable 1, the next pmax positions for
        % variable 2, etc
ordersinV=zeros(1,K);
% Initially the mean square error (MSE) of the DR model is equal to the
% variance of the response variable, and the BIC gets the largest value
MSEval=DRfitmse(xM,responseindex,ordersinV,indexV);
BICold=(n-max(ordersinV))*log(MSEval)+sum(ordersinV)*log(n-max(ordersinV));
ingameV=[1:K]; % Keep track of the variables for which the search has not
% yet reached the maximum order. If empty then terminate
ningame=K;
terminateflag=0;
incrisor=1; % The incremental step, starts with one and if none of the
% K candidate lags is selected, it becomes two etc.
while (~terminateflag && ningame~=0);
    pmaxreach=find(ordersinV>=pmax);
    ingameV=setdiff(ingameV,pmaxreach);
    ningame=length(ingameV);
    if ningame~=0
        BICnowV = NaN(ningame,1);
        MSEnowV = NaN(ningame,1);
        for iK=1:ningame
            ordtempV=ordersinV;
            ordtempV(ingameV(iK))=ordtempV(ingameV(iK))+incrisor;
            overpmaxV=find(ordtempV>pmax);
            if ~isempty(overpmaxV)
                ordtempV(overpmaxV)=pmax;
            end
            if length(overpmaxV)==K
                terminateflag=1;
            end
            tempindexV = indexV;
            tempindexV((ingameV(iK)-1)*pmax+ordtempV(ingameV(iK)))=1;
            MSEnowV(iK)=DRfitmse(xM,responseindex,ordtempV,tempindexV);
            BICnowV(iK)=(n-max(ordtempV))*log(MSEnowV(iK))+...
                sum(tempindexV)*log(n-max(ordtempV));
        end
        [BICnew,iBICnew]=min(BICnowV);
        invarindex=ingameV(BICnowV==BICnew);
        if (BICold<=BICnew)
            incrisor=incrisor+1;
            if incrisor>pmax-min(ordersinV)
                terminateflag=1;
            end
        else
            indexV((invarindex-1)*pmax+ordersinV(invarindex)+incrisor)=1;
            ordersinV(invarindex)=ordersinV(invarindex)+incrisor;
            BICold=BICnew;
            incrisor=1;
            MSEval = MSEnowV(iBICnew);
        end
    else
        terminateflag=1;
    end % if ningame==0
end % while terminate
% To facilitate further computations change the order of components in the
% index vector to first all variables at lag one, then all variables at lag
% two etc.
indexV=reshape(indexV,pmax,K);
indexV=indexV';
indexV=reshape(indexV,1,K*pmax);
maxorder=max(ordersinV);
end






