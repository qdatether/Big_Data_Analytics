function [RCGCIM,pRCGCIM]=mBTSCGCImatrix(xM,pmax,maketest)
% [RCGCIM,pRCGCIM] = mBTSCGCImatrix(xM,pmax,maketest)
% mBTSCGCImatrix computes the restricted conditional Granger Causality
% index (RCGCI) with the modified BTS method (mBTS) for all variable pairs
% (Xi,Xj) in the presence of the other obvserved variables. The $K$
% time series of the $K$ observed variables are given in the matrix 'xM'
% columnwise, and the maximum order for which mBTS searches for lagged 
% variables is given as 'pmax'. If selected ('maketest'=1), parametric 
% tests are also performed for each RCGCI. The parametric test is the 
% Fisher test with degrees of freedom for the Fisher distribution 
% determined by the mBTS regression model for each response variable Xj. 
% See also, Econometric Analysis, Greene, 7th Edition, Sec 5.5.2.
% For each response variable the function mBTSCGCI.m is called to run
% the mBTS and compute the RCGCI from each of the K-1 variables to the
% response variable.
% INPUTS
% - xM        : time series matrix (columns -> variables, rows ->
%               observations)
% - pmax      : maximum order for mBTS
% - maketest  : If 1 make parametric tests and give out the p-values
% OUTPUTS
% - RCGCIM    : the matrix of size KxK of the RCGCI values. Note that the
%               cell (i,j) corresponds to RCGCI(Xi->Xj).
% - pRCGCIM   : the matrix of size KxK of p-values of the parametric
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

nvar=size(xM,2);
RCGCIM=NaN(nvar,nvar);
pRCGCIM=NaN(nvar,nvar);
for j=1:nvar
    [RCGCIV,pRCGCIV] = mBTSCGCI(xM,j,pmax,maketest);
    RCGCIM(:,j)=RCGCIV;
    if maketest
        pRCGCIM(:,j)=pRCGCIV;
    end
end
