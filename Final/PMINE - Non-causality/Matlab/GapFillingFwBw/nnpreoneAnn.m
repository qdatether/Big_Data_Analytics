function xpre = nnpreoneAnn(xM,yV,tarV,nnei,q,isoutofsample)
% xpre = nnpreoneAnn(xM,yV,tarV,m,tau,nnei,q,isoutofsample)
% The function nnpreoneAnn takes in matrix of data points 'xM' and their 
% mappings 'yV', a query (target) point 'tarV' and the number of neighbors 
% 'nnei'. It gives to the output the predicted value for 'tarV' using the 
% 'nnei' nearest neighbor prediction. The type of prediction is specified
% by the input 'q' as follows:
% - q=0 : the average of the 'nnei' mappings. This is the default.
% - q=m, where m is the dimension for the points: the ordinary least square
%        solution for the linear model restricted to the 'nnei' neighbors.
% - 0<q<m : the principal component regression is applied to the 'nnei'
%        neighbors and the 'q' principal components are further used to
%        find the solution for the linear model.
% If the target point 'tarV' belongs to 'xM' then the input parameter
% 'isoutofsample' should be set to 0, if not it should be any other number. 
% Default is 1.

if nargin == 5
    isoutofsample = 1;
elseif nargin == 4
    isoutofsample = 1;
    q = 0;
elseif nargin == 3
    isoutofsample = 1;
    q = 0;
    nnei = 1;
end
if isempty(nnei), nnei=1; end
if isempty(q), q=0; end
if isempty(isoutofsample), isoutofsample=0; end

if isoutofsample
    neiindV = annquery(xM',tarV',nnei);
else
    neiindV = annquery(xM',tarV',nnei+1);
    neiindV = neiindV(2:nnei+1);
end
neiM = xM(neiindV,:);
yyV = yV(neiindV);
if q==0 || nnei==1
    xpre = mean(yyV);
else
    mneiV = mean(neiM);
    my = mean(yyV);
    zM = neiM - ones(nnei,1)*mneiV;
    [Ux, Sx, Vx] = svd(zM, 0);
    tmpM = Vx(:,1:q) * inv(Sx(1:q,1:q)) * Ux(:,1:q)';
    lsbV = tmpM * (yV - my);
    xpre = my + (tarV-mneiV) * lsbV;
end
