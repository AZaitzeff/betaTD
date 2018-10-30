function runMStd(filename,filesave,fid,numpar,betathres,Ks,dx,dy,dt,step,num)
addpath('../anglelib/')
if nargin<6
    Ks=[5,5];
end
if nargin<8
    dx=1/100;
    dy=1/100;
end
if nargin<9
    dt=.02;
end
if nargin<10
    step=1;
end
if nargin<11
    num=1;
end
parpool(numpar)

EBSDtemp=load(['../data/' filename 'EBSD.mat']);
addpath('../anglelib/')
EBSD=EBSDtemp.EBSD(rows,cols,:);
CI=EBSDtemp.CI(rows,cols);
%betas=EBSDtemp.betas(rows,cols);
reestbeta=1;
tic;
MStd(EBSD,CI,fid,Ks,filesave,dx,dy,dt,step,num,reestbeta);
toc;

poolobj = gcp('nocreate');
delete(poolobj);