function xnewM = fillgapsnnfbAnn(xM,noteM,tbeforeTMS,tafterTMS,taus,tbeforerec,tafterrec,memb,nnei)
% xnewM = fillgapsnnfbAnn(xM,noteM,tbeforeTMS,tafterTMS,taus,tbeforerec,tafterrec,memb,nnei)
% Same as fillgapsnnfb but using ANN for nearest neighbor search 
% INPUT
% - xM          : the EEG data, size N x m, N samples, m channels
% - noteM       : k x 2 matrix of annotations, k annotations, first column
%                 annotation code, second column second it appears.
% - tbeforeTMS  : start of corrupted signal with respect to TMS mark, in
%                 seconds.
% - tafterTMS   : end of corrupted signal with respect to TMS mark, in
%                 seconds.
% - taus        : the sampling time.
% - tbeforerec  : time before the TMS mark to which to search for neighbors
%                 and replace the corrupted signal, in seconds.
% - tafterrec   : time after the TMS mark to which to search for neighbors
%                 and replace the corrupted signal, in seconds.
% - memb        : embedding dimension to use for forming neighboring points
% - nnei        : number of neighbors to be used in smoothing
% OUTPUT
% - xnewM       : the smoothed EEG data matrix.

codeTMS = 1;  % the start of each TMS is given by this code
r = 1000;  % to deviate the range of the data and get the minimum distance
f = 1.2;  % factor to increase the distance if not enough neighbors are found
tau = 1; % Assume the delay is one.

[n,m]=size(xM);
iV = find(noteM(:,1)==codeTMS);
TMSnoteM = noteM(iV,:);
tV = [1:n]'*taus;
xnewM = xM;
sbeforeTMS = round(tbeforeTMS/taus);
safterTMS = round(tafterTMS/taus);
sbeforerec = round(tbeforerec/taus);
safterrec = round(tafterrec/taus);
fprintf('marks to replace segment=%d. Runinng... ',length(iV));
Tmax = safterTMS+sbeforeTMS+1;
for i=1:length(iV)
    if mod(i,50)==1
        fprintf('\n %d',i);
    else
        fprintf('.');
    end    
    [tmin,imin]=min(abs(tV-TMSnoteM(i,2)));
    for j=1:m
        x1V = xnewM([imin-sbeforerec:imin-sbeforeTMS-1],j);
        x2V = xnewM([imin+safterrec:-1:imin+safterTMS+1],j);
        % Forward prediction based on x1V
        nx1 = length(x1V);
        % State space reconstruction of the the training set
        nvec1 = nx1-1 - (memb-1)*tau;   % The length of the reconstructed set
        x1M = zeros(nvec1,memb);
        for iemb=1:memb
           x1M(:,memb-iemb+1) = x1V(1+(iemb-1)*tau:nvec1+(iemb-1)*tau);
        end
        y1V = x1V(2+(memb-1)*tau:nx1); % The one-step ahead mappings
        pre1V = NaN*ones(Tmax,1);
        winnowV = x1V(nx1-(memb-1)*tau:nx1);
        for T = 1:Tmax
            % Calls the function that makes one step prediction
            tarV = winnowV((memb-1)*tau+1:-tau:1)';    
            pre1V(T) = nnpreoneAnn(x1M,y1V,tarV,nnei,0,1);
            winnowV = [winnowV(2:end);pre1V(T)];
        end
        % Backward prediction based on x1V
        nx2 = length(x2V);
        % State space reconstruction of the the training set
        nvec2 = nx2-1 - (memb-1)*tau;   % The length of the reconstructed set
        x2M = zeros(nvec2,memb);
        for iemb=1:memb
           x2M(:,memb-iemb+1) = x2V(1+(iemb-1)*tau:nvec2+(iemb-1)*tau);
        end
        y2V = x2V(2+(memb-1)*tau:nx2); % The one-step ahead mappings
        pre2V = NaN*ones(Tmax,1);
        winnowV = x2V(nx2-(memb-1)*tau:nx2);
        for T = 1:Tmax
            % Calls the function that makes one step prediction
            tarV = winnowV((memb-1)*tau+1:-tau:1)';    
            pre2V(T) = nnpreoneAnn(x2M,y2V,tarV,nnei,0,1);
            winnowV = [winnowV(2:end);pre2V(T)];
        end
        % Give complementary weights according to their distance from start       
        xnewM(imin-sbeforeTMS:imin+safterTMS,j)=([Tmax-1:-1:0]'.*pre1V+[0:Tmax-1]'.*flipud(pre2V))/(Tmax-1);
    end
end
fprintf('\n');


