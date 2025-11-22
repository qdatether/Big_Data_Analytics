function [RCGCIV,pRCGCIV] = mBTSCGCI(xM,responseindex,pmax,maketest)
% [RCGCIV,pRCGCIV] = mBTSCGCI(xM,responseindex,pmax,maketest)
% mBTSCGCI computes the restricted conditional Granger Causality
% index (RCGCI) with the modified BTS method (mBTS) for all variable pairs
% (Xi,Xj) in the presence of the other obvserved variables, where Xj is
% the given response variable indicated by 'responseindex'. The $K$
% time series of the $K$ observed variables are given in the in 'xM'
% columnwise and 'responseindex' gives the column order of the response
% variable. If selected, parametric tests are also performed for each
% of the K-1 RCGCI. The parametric test is the Fisher test with degrees of
% freedom for the Fisher distribution determined by the mBTS regression
% model for each response variable Xj. See also, Econometric Analysis,
% Greene, 7th Edition, Sec 5.5.2.
% INPUTS
% - xM        : time series matrix of K columns (columns -> variables, rows
%               -> observations)
% - responseindex: the index for the response, an integer in {1,...,K}
% - pmax      : maximum order for mBTS
% - maketest  : If 1 make parametric tests and give out the p-values
% OUTPUTS
% - RCGCIV    : vector of size Kx1 of the RCGCI values. Note that the
%               component at position 'responseindex' is NaN.
% - pRCGCIM   : the vector of size Kx1 of p-values of the parametric
%               significance test for the corresponding RCGCI values.
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

if nargin==3
    maketest = 0;
end
[n,K]=size(xM);
RCGCIV = zeros(K,1);
pRCGCIV = NaN(K,1);
[indexV,maxorder]=mBTS(xM,responseindex,pmax);
if maxorder>0
    % If maxorder==0 the unrestricted model is empty, so the response has
    % no auto- and cross- correlation, thus CGCI=0.
    RCGCIV(responseindex) = NaN;
    for d=1:K
        xM(:,d)=xM(:,d)-mean(xM(:,d));
    end
    tmpindexV = indexV(1:maxorder*K);
    % In order to feed it to DRfitnrmse the sequence of components should
    % change to first all lags of first variable, then all lags of second
    % variable etc
    unrestrictedindexV=reshape(tmpindexV,K,maxorder);
    unrestrictedindexV=unrestrictedindexV';
    unrestrictedindexV=reshape(unrestrictedindexV,1,K*maxorder);
    MSEunrestricted=DRfitmse(xM,responseindex,maxorder*ones(1,K),...
        unrestrictedindexV);
    allactive = sum(unrestrictedindexV);
    drivingV = setdiff([1:K],responseindex);
    for iK=drivingV
        % Dynamic regression on Y from Y-past and Z-past excluding X-past
        restrictedindexV = unrestrictedindexV;
        restrictedindexV((iK-1)*maxorder+1:iK*maxorder) = zeros(1,maxorder);
        MSErestricted=DRfitmse(xM,responseindex,maxorder*ones(1,K),...
            restrictedindexV);
        RCGCIV(iK) = log(MSErestricted/MSEunrestricted);
        if maketest
            % Compute the p-value of the parametric F-test for the significance of
            % RCGCI
            if RCGCIV(iK)==0
                % This is the case no components of the driving variable are in
                % the model of mBTS
                pRCGCIV(iK) = 1;
            else
                parnumerator=sum(unrestrictedindexV((iK-1)*maxorder+1:iK*maxorder));
                pardenominator = n-maxorder-allactive;
                if pardenominator>0
                    xfstat = ((MSErestricted-MSEunrestricted)/parnumerator)...
                        /(MSEunrestricted/pardenominator);
                    pRCGCIV(iK) = 1 - fcdf(xfstat,parnumerator,pardenominator);
                else
                    pRCGCIV(iK) = 1;
                    fprintf('WARNING: Cannot make the parametric test for driving variable No %d.\n',iK);
                    fprintf('         Number of explanatory variables >= number of effective data points.');
                end
            end
        end % if maketest
    end % for iK
end % if maxorder

        
        
