function MSE=DRfitmse(xM,responseindex,ordersV,indexV)
% DRfitmse makes predictions with the given dynamic regression model,
% determined by the lag indexes for all K variables in 'indexV'. The data
% matrix 'xM' contains all K variables columnwise, the response variables
% is given by 'responseindex', and the maximum order of each of the K 
% variables is given in 'ordersV'. 
% INPUTS
% - xM          : the matrix of $K$ time series (variables in columns)
% - responseindex : the index of the response variable in {1,...,K}
% - ordersV     : vector of size 1xK of the maximum order for each of the K
%                 variables.
% - indexV      : the vector of size 1 x K*pmax of zeros and ones
%                 e.g. if the component in position 2*K+3 is one, the
%                 second variable, lag 3, X2(t-3), is selected.
% OUTPUTS
% - MSE         : the mean square error (MSE) of the fit of the dynamic
%                 regression model.
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

n = size(xM,1);
xlagM = multilagmatrix(xM,responseindex,ordersV,indexV);
An=inv(xlagM(:,2:end)'*xlagM(:,2:end))*xlagM(:,2:end)'*xlagM(:,1);
preV = xlagM(:,2:end)*An;
resV = xM(max(ordersV)+1:end,responseindex)-preV;
MSE = (resV'*resV/length(resV));
