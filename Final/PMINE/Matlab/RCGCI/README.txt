4/8/2015, Thessaloniki, Elsa Siggiridou and Dimitris Kugiumtzis
The folder contains the Matlab codes to run restricted conditional Granger causality index (RCGCI) based on the dynamic regression model derived by the modified backward-in-time selection (mBTS) algorithm. 
The folder contains the following files:
- Example.m : a script on an exemplary multivariate time series from a known multivariate stochastic process on K variables.
- mBTSCGCImatrix.m : the main function called in the script, giving out the matrix (size K x K) of RCGCI as well as the matrix of the corresponding p-values from the parametric significance test.
- mBTSCGCI.m : the function called by 'mBTSCGCImatrix' that computes the RCGCI an the corresponding p-values for a given response variable (and the K-1 driving variables).
- mBTS.m : the function called in 'mBTSCGCI' to run mBTS for the given response variable.
- DRfitmse.m : the function called in 'mBTS' and also in 'mBTSCGCI' to fit the given dynamic regression model and give out the mean square error (MSE).
- multilagmatrix.m : called in 'DRfitmse' to build the matrix of explanatory variables.
- SysSchelter06bm1.m : the function generating multivariate time series from a known multivariate stochastic process on K variables. 