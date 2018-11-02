function runMStd(filename,filesave,fid,numpar,betathres,Ks,dx,dy,dt,step,num,reestbeta,mexed)
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
    reestbeta=0;
end
if nargin<13
    mexed=0;
end

parpool(numpar)

EBSDtemp=load(['../data/' filename 'EBSD.mat']);
addpath('../anglelib/')
EBSD=EBSDtemp.EBSD;
CI=EBSDtemp.CI;
[m,n,z]=size(EBSD);
%betas=EBSDtemp.betas(rows,cols);

tic;
parfor i=1:num
MStd(EBSD,CI,fid,Ks,filesave,dx,dy,dt,step,i,mexed);
end
toc;
poolobj = gcp('nocreate');
delete(poolobj);

enevec=zeros(1,num);
for i=1:num
    vars=load(['results/' filesave num2str(fid) num2str(i)]);
    enevec(i)=vars.energy;
end
[~,I]=min(enevec);
vars=load(['results/' filesave num2str(fid) num2str(I)]);


mapall=vars.mapall;
dict=vars.dict;
energy=vars.energy;
if reestbeta
    [dict,~]=estimatebetasfast(EBSD,CI,mapall,24,mexed);
end
betaEBSD=zeros(m*n,4);
K = max(mapall(:));
alpha=reshape(EBSD, [m*n,z]);
alpha=E313toq(alpha);
for z=1:K
    mu=dict{z};
    zfun=find(z==mapall);
    [best,~]=alpbbeta(alpha(zfun,:),mu);
    betaEBSD(zfun,:)=best;

end
betaEBSD=qtoE313(betaEBSD);
betaEBSD=reshape(betaEBSD,[m,n,3]);
save(['results/' filesave num2str(fid)],'mapall','betaEBSD','dict','energy');

for i=1:num
    delete(['results/' filesave num2str(fid) num2str(i) '.mat']);
end

