% Fit an ARMA model. 
% 1. Generation of a time series from a given ARMA process, or load a given
% time series
% 2. Autocorrelation and partial autocorrelation 
% 3. Fit of a ARMA model.
% 4. Prediction error statistic for the ARMA model.
% 5. Multi-step predictions from a given target point with the ARMA model.
clear all
rooV = [0.4+0.3*i 0.4-0.3*i 0.6 0.8]';
% rooV = [0.8 0.9 0.6 0.8]';
phi0= 0;
thetaV = []; % [-5 0.4];
n = 100;
sdnoise = 1;
maxtau = 10;
pmax = 6;
qmax = 3;
Tmax = 10;
proptest = 0.3;
alpha = 0.05;

tmpV = poly(rooV);
phiV = [phi0; -tmpV(2:length(rooV)+1)'];
pgen = length(phiV)-1;
qgen = length(thetaV);
% 1. Generation of a time series from a given ARMA process
xV = generateARMAts(phiV,thetaV,n,sdnoise);
stext = sprintf('ARMA(%d,%d)',pgen,qgen);
% stext = 'RRMA1';
% xV = load(strcat(stext,".dat"));
n = length(xV);

figure(1)
clf
plot(xV,'.-')
hold on
xlabel('t')
ylabel('x(t)')
title(sprintf('%s, time history',stext))

% 2. Autocorrelation and partial autocorrelation 
% 2a. Autocorrelation
[acM] = autocorrelation(xV, maxtau);
zalpha = norminv(1-alpha/2);
autlim = zalpha/sqrt(n);
figure(2)
clf
hold on
for ii=1:maxtau
    plot(acM(ii+1,1)*[1 1],[0 acM(ii+1,2)],'b','linewidth',1.5)
end
plot([0 maxtau+1],[0 0],'k','linewidth',1.5)
plot([0 maxtau+1],autlim*[1 1],'--c','linewidth',1.5)
plot([0 maxtau+1],-autlim*[1 1],'--c','linewidth',1.5)
xlabel('\tau')
ylabel('r(\tau)')
title(sprintf('%s, autocorrelation',stext))

% 2b. Partial autocorrelation
display = 1;
pacfV = parautocor(xV,maxtau);
figure(3)
clf
hold on
for ii=1:maxtau
    plot(acM(ii+1,1)*[1 1],[0 pacfV(ii)],'b','linewidth',1.5)
end
plot([0 maxtau+1],[0 0],'k','linewidth',1.5)
plot([0 maxtau+1],autlim*[1 1],'--c','linewidth',1.5)
plot([0 maxtau+1],-autlim*[1 1],'--c','linewidth',1.5)
xlabel('\tau')
ylabel('\phi_{\tau,\tau}')
title(sprintf('%s, partial autocorrelation',stext))

aicM = NaN(pmax+1,qmax+1);
for p=0:pmax
    for q=0:qmax
        [nrmseV,phiallV,thetaallV,SDz,aicM(p+1,q+1),fpeS]=fitARMA(xV,p,q,1);
    end
end
figure(4)
clf
plot([0:pmax],aicM(:,1),'.-k','Markersize',12,'linewidth',1.5)
hold on
if qmax>=1
    plot([0:pmax],aicM(:,2),'.-b','Markersize',12,'linewidth',1.5)
end
if qmax>=2
    plot([0:pmax],aicM(:,3),'.-c','Markersize',12,'linewidth',1.5)
end
if qmax>=3
    plot([0:pmax],aicM(:,4),'.-g','Markersize',12,'linewidth',1.5)
end
if qmax>=4
    plot([0:pmax],aicM(:,5),'.-r','Markersize',12,'linewidth',1.5)
end
if qmax>=5
    plot([0:pmax],aicM(:,6),'.-y','Markersize',12,'linewidth',1.5)
end
xlabel('p')
ylabel('AIC(p,q)')
title(sprintf('%s, AIC',stext))
if qmax == 1
    legend('q=0','q=1')
end
if qmax == 2
    legend('q=0','q=1','q=2')
end
if qmax == 3
    legend('q=0','q=1','q=2','q=3')
end
if qmax == 4
    legend('q=0','q=1','q=2','q=3','q=4')
end
if qmax == 5
    legend('q=0','q=1','q=2','q=3','q=4','q=5')
end

% 3. Fit of an ARMA model.
p = input('Give the order p of the AR part >'); 
q = input('Give the order q of the MA part >'); 
[nrmseV,phiallV,thetaallV,SDz,aicS,fpeS]=fitARMA(xV,p,q,1);
fprintf('===== ARMA model ===== \n');
fprintf('Estimated coefficients of phi(B):\n');
disp(phiallV')
fprintf('Estimated coefficients of theta(B):\n');
disp(thetaallV')
fprintf('SD of noise: %f \n',SDz);
fprintf('AIC: %f \n',aicS);
fprintf('FPE: %f \n',fpeS);
fprintf('NRMSE: %f \n',nrmseV);

% 4. Prediction error statistic for the ARMA model.
nlast = round(proptest*n);
tittxt = sprintf('ARMA(%d,%d), %%test=%1.2f, prediction error',p,q,proptest);
figure(5);
clf
[nrmseV,preM] = predictARMAnrmse(xV,p,q,Tmax,nlast,'example');
figure(6);
clf
plot([n-nlast+1:n]',xV(n-nlast+1:n),'.-')
hold on
plot([n-nlast+1:n]',preM(:,1),'.-r')
if Tmax>1
    plot([n-nlast+1:n]',preM(:,2),'.-c')
	if Tmax>2
        plot([n-nlast+1:n]',preM(:,3),'.-k')
    end
end
switch Tmax
    case 1
        legend('true','T=1','Location','Best')
    case 2
        legend('true','T=1','T=2','Location','Best')
    otherwise
        legend('true','T=1','T=2','T=3','Location','Best')
end
% 5. Multi-step predictions from a given target point with the ARMA model.
n1 = n-Tmax;
figure(7);
clf
[preV] = predictARMAmultistep(xV,n1,p,q,Tmax,'example');
