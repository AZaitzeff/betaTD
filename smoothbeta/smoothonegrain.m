addpath('../anglelib/');

load('testEthan.mat');


[m,n,~]=size(EBSD);
smallflatq=E313toq(reshape(EBSD,m*n,3)');
[mu, kap] = estimatebeta(m*n, CI(:), smallflatq, beta(:), [], [], 400);
propbeta = RMatOfQuat(mu);
Rmat = bestparentmatrix(m,n,smallflatq,beta,propbeta,1);

%your code here
smoothed = smoothBetaInterp(propbeta, Rmat);
