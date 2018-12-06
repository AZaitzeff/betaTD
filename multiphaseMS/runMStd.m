function runMStd(filename,filesave,fid,numpar,betathres,Ks,dx,dy,dt,step,num,reestbeta,clean,subsample,enec)
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
    clean=-1;
end
if nargin<14
    enec=0;
end

EBSDtemp=load(['../data/' filename 'EBSD.mat']);
addpath('../anglelib/')
EBSD=EBSDtemp.EBSD;
CI=EBSDtemp.CI;
[M,N]=size(CI);
beta=logical(EBSDtemp.betas);
clean=min(clean,.9);
K=prod(Ks);
%codegenzaitzeff(M,N,K);
if clean>0
    [EBSD,CI]=cleanEBSDdata(EBSD,CI,clean);
end
if subsample>1
    EBSD=EBSD(1:subsample:end,1:subsample:end,:);
    CI=CI(1:subsample:end,1:subsample:end);
    beta=beta(1:subsample:end,1:subsample:end);
end

[m,n,z]=size(EBSD);
%betas=EBSDtemp.betas(rows,cols);

tic;
if numpar>1
    parpool(numpar)
    parfor i=1:num
        MStd(EBSD,CI,beta,fid,Ks,filesave,dx,dy,dt,i,enec);
    end
    poolobj = gcp('nocreate');
    delete(poolobj);
    
else
    for i=1:num
        MStd(EBSD,CI,beta,fid,Ks,filesave,dx,dy,dt,i,enec);
    end
end
toc;


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
for k=1:K
    mu=dict(k,:);
    zfun=find(k==mapall);
    [best,~]=alpbbeta(alpha(zfun,:),mu);
    betaEBSD(zfun,:)=best;

end
betaEBSD=qtoE313(betaEBSD);
betaEBSD=reshape(betaEBSD,[m,n,3]);
save(['results/' filesave num2str(fid)],'mapall','betaEBSD','dict','energy');

for i=1:num
    delete(['results/' filesave num2str(fid) num2str(i) '.mat']);
end

