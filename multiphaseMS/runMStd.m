function runMStd(filename,filesave,fid,numpar,betathres,Ks,dx,dy,dt,step,num,mexed)
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
if nargin<12
    mexed=0;
end

parpool(numpar)

EBSDtemp=load(['../data/' filename 'EBSD.mat']);
addpath('../anglelib/')
EBSD=EBSDtemp.EBSD;
CI=EBSDtemp.CI;
%betas=EBSDtemp.betas(rows,cols);
reestbeta=0;
tic;
MStd(EBSD,CI,fid,Ks,filesave,dx,dy,dt,step,num,reestbeta,mexed);
toc;

poolobj = gcp('nocreate');
delete(poolobj);