% The script runs the RCGCI with mBTS on time series from a VAR-process
% of order 4 in 5 variables.
N=1000; % time series length
maxord=6; % maximum order for mBTS
maketest=1; % If 1 make parametric tests and give out the p-values

rand('state',1);
randn('state',1);
xM=SysSchelter06bm1(N); %multivariate time-series
[RCGCIM,pRCGCIM] = mBTSCGCImatrix(xM,maxord,maketest)