addpath('../anglelib/')
load(['../data/' 'AFone' 'EBSD'])
 [M,N]=size(IQ);
% midm=ceil(M/2);
% midn=ceil(N/2);
% filler=200;
% lm=max(midm-filler,1);
% um=min(midm+filler,M);
% ln=max(midn-filler,1);
% un=min(midn+filler,N);
% suby=lm:um;
% subx=ln:un;
% M=um-lm+1;
% N=un-ln+1;
% EBSD=EBSD(suby,subx,:);
% IQ=IQ(suby,subx);
% CI=CIfunc(CI(suby,subx));
% beta=logical(betas(suby,subx));
 codegenpilchank(M,N)
 CI=CIfunc(CI);
  beta=logical(betas);


tic;[newmap,newbetadict,map,betadict,K,filled]=fullpilchank_mex(EBSD,beta,IQ,CI,5,5,5,5,5);toc;
[EBSD,valued]=map2bungecoord(map,betadict,filled,K);

%addpath('../../../MATLAB/mtex-5.1.1/')
%startup

%[replacemap,betadict,kappa,numgrains]=Germain(neigharray,sizeneigh,dict,betamean,IQmean,K,5,5);
%newmap=changemapleavezeros(map,replacemap);