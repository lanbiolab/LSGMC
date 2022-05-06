clc
clear
addpath('dataset');addpath('other function');
load('COIL20MV');
lambda = 10;%[1e-5 1e-4 1e-3 0.01 0.1 1 10 100]
ppp=0.5;%[0.01 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9]
%parameter  lambda、ppp
%MSRC-v1: 0.1 0.1      yaleB,handwritten,COIL20MV: 10 0.5       flower17: 1e-3 0.7
%proteinFold: 10 0.9   3sources: 100 0.9                        Caltech101-7: 10 0.6
%% 
for i=1:size(X,2)
    X{i} = X{i}./(repmat(sqrt(sum(X{i}.^2,1)),size(X{i},1),1)+eps);
end
nCluster = length(unique(Y));
mu=1e-2;
rho=2.8;

tic
Zn = LSGMC(X,lambda,mu,rho,ppp);
M=retain(Zn);
W=postprocessor(M);
label = new_spectral_clustering(W,nCluster);
[acc,nmi,F,precision,AR,~,~] = AllMeasure(label,Y);
disp(['acc   ' num2str(acc*100) '%     ' 'nmi   ' num2str(nmi*100) '%'  '     time     '   num2str(toc)])
