addpath('../betaTD')
addpath('../anglelib/')
numsub=400;
dx=1/100;
dy=dx;
dtstop=2^(-10);
between=2;
nt=8;
DT=2^-4;
load('testmultigrain.mat');
CI=CIfunc(CI);




[M,N]=size(mapall);



 vectorType1 = coder.typeof(1, [M N], [false false]);
 vectorType3 = coder.typeof(1, [4 M*N], [false false]);
 vectorType1l = coder.typeof(1==1, [M N], [false false]);
 vectorType2 = coder.typeof(1, [M N 3], [false false false]);
 
 codegen EBSDMStdsmooth -args {vectorType1,vectorType2,vectorType1,vectorType1l,vectorType3,1,1,1,1,1,1,1,1,1,1}

%%
fid=75;
dtsm=1/1000;
ntsm=100;
[mapall2,smoothEBSD2,energy,flag]=EBSDMStdsmooth_mex(mapall,EBSD,CI,beta,smoothEBSD,K,fid,DT,dx,dy,dtstop,nt,between,ntsm,dtsm);
plotsmoothwbnds(mapall,smoothEBSD)
figure
plotsmoothwbnds(mapall2,smoothEBSD2)