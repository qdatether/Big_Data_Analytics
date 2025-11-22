function [PTERVxy,PTERVyx] = PTERV(xV,yV,zM,T,mx,taux,my,tauy,mz,tauz)
% [PTERVxy,PTERVyx] = PTERV(xV,yV,zM,T,mx,taux,my,tauy,mz,tauz)
% PTERV estimates the partial transfer entropy on rank vectors (PTERV). 
% PTERV indicates the direct influence of time series xV on time series yV 
% accounting for the influence on yV from other time series contained in
% zM. Then the same measure is computed for the direct driving of xV from 
% yV and the resulting values PTERVxy and PTERVyx, respectively, are given
% to the output. 
% PTERV is the extension of TERV for multivariate time series, 
% D. Kugiumtzis "Transfer Entropy on Rank Vectors", Journal of Nonlinear
% Systems and Applications, Vol 3, No 2, pp 73-81, 2012 [arxiv:1007.0394] 
% Note that the embeddings of X and Y are the same for both directions.
% The definition of TERV is the same as for partial transfer entropy (PTE) 
% PTERVxy = PTERV(X->Y)= 
% -H([y(t+T)|xk(t),yl(t),zm(t)]) + H([y(t+T)|yl(t),zm(t)]) = 
%  -H([y(t+T),xk(t),yl(t),zm(t)]) + H([xk(t),yl(t),zm(t)])
%  +H([y(t+T),yl(t),zm(t)]) - H([yl(t),zm(t)])
% 
% where xk(t) = ranks of [x(t),x(t-taux),...,x(t-(mx-1)taux)]'
%       yl(t) = ranks of [y(t),y(t-tauy),...,y(t-(my-1)tauy)]'
%       zm(t) = ranks of [z(t),z(t-tauz),...,y(t-(mz-1)tauz)]'
% The embedding of more than one conditioning variables is done as for
% zm(t).
% INPUTS
% - xV   : time series 1, column vector
% - yV   : time series 2, column vector
% - zM   : the other time series, columnwise 
% - T    : T steps ahead, note that if T>1 the whole future vector of
%          length T is considered.
% - mx   : embedding dimension for xV
% - taux : lag for xV
% - my   : embedding dimension for yV 
% - tauy : lag for yV
% - mz   : embedding dimension for time series in zM, the same embedding
%          dimension is used for all time series in zM. 
% - tauz : lag for time series in zM, as for mz
% OUTPUTS
% - PTERVxy: partial transfer entropy from X to Y 
% - PTERVyx: partial transfer entropy from Y to X, both for the same
%            embeddings of X and Y 
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
% Reference : D. Kugiumtzis, "Partial Transfer Entropy on Rank Vectors", 
% The European Physical Journal Special Topics, Vol 222, No 2, pp 401-420, 
% 2013
% Link      : http://users.auth.gr/dkugiu/
%=========================================================================

N = size(xV,1);
if size(yV,1)~=N || size(zM,1)~=N || size(yV,2)~=1 || size(xV,2)~=1 
    error('Not proper data matrix sizes.');
end
nz = size(zM,2);
M = max([(mx-1)*taux (my-1)*tauy (mz-1)*tauz]);
N1 = N-M-T;
% The rank joint vector of a) the future vector of the response X, 
% [x(t+1), x(t+2),...,x(t+T)], and b) the embedding vector of X,
% [x(t), x(t-taux), ..., x(t-(mx-1)taux)].
xfM = NaN*ones(N1,mx+T);
yfM = NaN*ones(N1,my+T);
for imx=1:mx
    xfM(:,imx) = xV(M+1-(imx-1)*taux:N-T-(imx-1)*taux);
end
for imy=1:my
    yfM(:,imy) = yV(M+1-(imy-1)*tauy:N-T-(imy-1)*tauy);
end
for iT=1:T
    xfM(:,mx+iT)=xV(M+1+iT:N-T+iT);
    yfM(:,my+iT)=yV(M+1+iT:N-T+iT);
end
[tmpM,rxfM]=sort(xfM,2); 
[tmpM,ryfM]=sort(yfM,2); 
[tmpM,rxfM]=sort(rxfM,2); 
[tmpM,ryfM]=sort(ryfM,2); 
% From the matrix of first mx columns compute the ranks, these are the
% ranks of the embedding vector of X alone.
[tmpM,rxM]=sort(rxfM(:,1:mx),2); % The ranks for the embedding vector of X
% [tmpM,rxM]=sort(rxM,2); 
[tmpM,ryM]=sort(ryfM(:,1:my),2);  % The ranks for the embedding vector of Y
% [tmpM,ryM]=sort(ryM,2);
rxfM = rxfM(:,mx+1:mx+T); % The ranks for the future vector of X
ryfM = ryfM(:,my+1:my+T); % The ranks for the future vector of X
% The matrix of words (at each row) for system Z
zzM = NaN*ones(N-M-T,nz*mz);
rzM = NaN*ones(N-M-T,nz*mz);
for inz=1:nz
    for imz=1:mz
        zzM(:,(inz-1)*mz+imz) = zM(M+1-(imz-1)*tauz:N-T-(imz-1)*tauz,inz);
    end
    [tmpM,rzM(:,(inz-1)*mz+1:(inz-1)*mz+mz)]=sort(zzM(:,(inz-1)*mz+1:(inz-1)*mz+mz),2);
end
% [tmpM,ryfM]=sort(ryfM,2); 
% % PTERVxy = PTERV(X->Y)= 
%  -H([ry(t+T),rxk(t),ryl(t),rzm(t)]) + H([rxk(t),ryl(t),rzm(t)])
%  +H([ry(t+T),ryl(t),rzm(t)]) - H([ryl(t),rzm(t)])
% H([rxk(t),ryl(t),rzm(t)]), where rx is the rank of x, ry is the rank of
% y, and rz is the rank of z 
H1 = symbolvectorentropy([ryM rxM rzM]);
% Direction X->Y
% H([y(t+T),xk(t),yl(t),zm(t)]) = H([ryl(t+T),rxk(t),rzm(t)]), where
% ryl(t+T) is the rank of vector 
% [y(t-(my-1)*tauy) y(t-(my-2)*tauy) ... y(t), y(t+T)]
H2xy = symbolvectorentropy([ryfM ryM rxM rzM]);
% H([y(t+b),yl(t),zm(t)])
H3xy = symbolvectorentropy([ryfM ryM rzM]);
% H([yl(t),zm(t)])
H4xy = symbolvectorentropy([ryM rzM]);
PTERVxy = - H2xy + H1 + H3xy - H4xy;
% Direction Y->X
H2yx = symbolvectorentropy([rxfM rxM ryM rzM]);
H3yx = symbolvectorentropy([rxfM rxM rzM]);
H4yx = symbolvectorentropy([rxM rzM]);
PTERVyx = - H2yx + H1 + H3yx - H4yx;

function H = symbolvectorentropy(xM,nsymbols)
% function H = symbolvectorentropy(xM)
% SYMBOLVECTORENTROPY computes the entropy of symbol vectors given in 'xM'.
% The replicates of each symbol vector in matrix 'xM' are counted to give
% the relative frequency, which is the estimate of the (non-zero)
% probability of the symbol vector. Then the entropy is computed from the
% probabilities of all different symbol vectors in 'xM'. 
% NOTE that 'xM' is supposed to contain natural numbers corresponding to
% symbols.
% INPUT
% - xM  : nxm matrix, each row is an m-dimensional array of natural numbers
% - nsymbols : if specified, it denotes the number of symbols in 'xM',
%              supposing that the symbols in 'xM' are 1,2,...,nsymbols. If 
%              not specified then the maximum natural number of 'xM' is 
%              found, and it is then supposed that all other components of 
%              xM are also natural numbers.
% OUTPUT
% - H   : the computed entropy.

if nargin==1
    nsymbols = 0;
end

[n,m] = size(xM);
% If not specified, compute the maximum natural number in 'xM'.
if nsymbols==0
    nsymbols = max(max(xM));
    if uint32(nsymbols)~=nsymbols
        error('The maximum symbol in the input matrix is not a natural number');
    end
end
% One-to-one transformation of each row of 'xM' to a number 
xV = xM*(nsymbols.^[0:m-1]');
clear xM
oxV = sort(xV);
doxV = diff(oxV);
uniqV = doxV ~= 0;
uniqV = [true; uniqV];
uniqV(n+1)=true;
ordindV = [1:n+1]';
ordindV = ordindV(uniqV);
pV = diff(ordindV);
pV = pV/n;  
H = -sum(pV.*log2(pV));
