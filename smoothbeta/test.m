addpath('../betaTD')
numsub=400;
dx=1/100;
dy=dx;
dtstop=2^(-10);
between=2;
nt=8;
DT=2^-4;
fid=75;
load('testmultigrain.mat');
CI=CIfunc(CI);

dtsm=1/10000;
ntsm=100;


[M,N]=size(mapall);
[~,~,z]=size(EBSD);
EBSDflat=E313toq(reshape(EBSD,[M*N,z])');
[xbdcor,ybdcor,sizebdcor,coords,sizecoords,minmaxrowcol]=  bndcoords(mapall,K);%Gets coordinates
k=2; %k goes from 1 to 14
csize=sizecoords(k);
indices=coords(k,1:csize);
[u,uin,wsm]=estso3smoothonegrain(M,N,indices,csize,EBSDflat,beta,CI,smoothEBSD,ntsm,dtsm,fid);
[im2, oM]=visso3(uin);
[im, oM]=visso3(u);
imagesc([im im2])

% vectorType1 = coder.typeof(1, [M N], [false false]);
% vectorType3 = coder.typeof(1, [4 M*N], [false false]);
% vectorType1l = coder.typeof(1==1, [M N], [false false]);
% vectorType2 = coder.typeof(1, [M N 3], [false false false]);
% 
% codegen EBSDMStdsmooth -args {vectorType1,vectorType2,vectorType1,vectorType1l,vectorType3,1,1,1,1,1,1,1,1,1,1}
% 
% [mapall2,smoothEBSD2,energy,flag]=EBSDMStdsmooth_mex(mapall,EBSD,CI,beta,smoothEBSD,K,fid,DT,dx,dy,dtstop,nt,between,ntsm,dtsm);