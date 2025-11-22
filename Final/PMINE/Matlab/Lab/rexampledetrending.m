% Detrending of a given time series assuming to have a trend component.
% Remove trend using moving average smoothing
clear all
n = 1000;
sdnoise = 10;
mux = 20;
maxtau = 100;
alpha = 0.05;
zalpha = norminv(1-alpha/2);
% 1. Generation of a time series with trend
xV = sdnoise * randn(n,1) + mux; % Gaussian iid time series
% xV = generateARMAts([0 0.5],[],n,1);
yV = addstochastictrend(xV); % add stochastic trend
% 2. Load a time series with trend
% stext = 'RR';
% yV = load(strcat(stext,".dat"));

%%%%%% Remove trend using a moving average filter
figure(1)
clf
plot(yV,'.-')
hold on
xlabel('t')
ylabel('y(t)')
title('time series with trend')

% Autocorrelation of the original time series 
acyM = autocorrelation(yV, maxtau);
autlim = zalpha/sqrt(n);
figure(2)
clf
hold on
for ii=1:maxtau
    plot(acyM(ii+1,1)*[1 1],[0 acyM(ii+1,2)],'b','linewidth',1.5)
end
plot([0 maxtau+1],[0 0],'k','linewidth',1.5)
plot([0 maxtau+1],autlim*[1 1],'--c','linewidth',1.5)
plot([0 maxtau+1],-autlim*[1 1],'--c','linewidth',1.5)
xlabel('\tau')
ylabel('r(\tau)')
title(sprintf('original time series, autocorrelation'))

maorder = input('Remove trend: Give an order for the moving average filter: > ');
mu1V = movingaveragesmooth2(yV,maorder);
figure(3)
clf
plot(yV,'.-')
hold on
plot(mu1V,'.-r')
xlabel('t')
ylabel('y(t)')
title('time series with trend')
legend('original',sprintf('MA(%d) smooth',maorder),'Location','Best')

x1V = yV - mu1V;
figure(4)
clf
plot(x1V,'.-')
xlabel('t')
ylabel('x(t)')
title(sprintf('detrended time series by MA(%d) smooth',maorder))

% Autocorrelation of the detrended time series by moving average filter 
acx1M = autocorrelation(x1V, maxtau);
figure(5)
clf
hold on
for ii=1:maxtau
    plot(acx1M(ii+1,1)*[1 1],[0 acx1M(ii+1,2)],'b','linewidth',1.5)
end
plot([0 maxtau+1],[0 0],'k','linewidth',1.5)
plot([0 maxtau+1],autlim*[1 1],'--c','linewidth',1.5)
plot([0 maxtau+1],-autlim*[1 1],'--c','linewidth',1.5)
xlabel('\tau')
ylabel('r(\tau)')
title(sprintf('detrended time series by MA(%d) smooth, autocorrelation',maorder))

