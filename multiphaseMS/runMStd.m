function runMStd(filename,filesave,numpar,num,gs)
if nargin<5
    gs=50;
end
addpath('../anglelib/')


if nargin<3
    numpar=10;
    num=10;
end


EBSDtemp=load(['../data/' filename 'EBSD.mat']);
addpath('../anglelib/')
dx=1/(2*gs);
dy=dx*EBSDtemp.scale;

[M,N]=size(EBSDtemp.CI);
midm=ceil(M/2);
midn=ceil(N/2);
lm=max(midm-250,1);
um=min(midm+250,M);
ln=max(midn-250,1);
un=min(midn+250,N);
suby=lm:um;
subx=ln:un;
M=um-lm+1;
N=un-ln+1;
EBSD=EBSDtemp.EBSD(suby,subx,:);
CI=EBSDtemp.CI(suby,subx);
beta=logical(EBSDtemp.betas(suby,subx));
codegenzaitzeff(M,N,1);

%betas=EBSDtemp.betas(rows,cols);
dt=2^-5;

nr=ceil(M/30*gs/50);
nc=ceil(N/30*gs/50);
[mapallp,dictp,kappap,~]=initializeEBSDfast_mex(EBSD,CI,beta,nr,nc);
K=nr*nc;
[~,~,~,coords,sizecoords,~]=  bndcoords(mapallp,K);
alpha=reshape(EBSD, [M*N,3]);
alpha=E313toq(alpha);
truebetaEBSD=zeros(M*N,4);
for k=1:K
    mu=dictp(k,:);
    indices=coords(k,1:sizecoords(k));
    mask=beta(indices);
    acoord=indices(~mask);
    [best,~]=alpbbeta(alpha(acoord,:),mu);
    truebetaEBSD(acoord,:)=best;
    bcoord=indices(mask);
    truebetaEBSD(bcoord,:)=alpha(bcoord,:);
end

nr=ceil(M/40*gs/50);
nc=ceil(N/40*gs/50);

fids=[50,100,150,200,250,300];
numfids=numel(fids);
runcheck=9;
totalcheck=runcheck*numfids;
if numpar>1
    parpool(numpar)
    
    parfor pari=1:totalcheck
        fidz=fids(ceil(pari/runcheck));
        MStd(EBSD,CI,beta,fidz,filesave,dt,dx,dy,nr,nc,mod(pari-1,runcheck)+1);
    end
    
    poolobj = gcp('nocreate');
    delete(poolobj);
    
else
    for i=1:totalcheck
        fid=fids(ceil(i/runcheck));
        MStd(EBSD,CI,beta,fid,filesave,dt,dx,dy,nr,nc,mod(i-1,runcheck)+1);
    end
end
total=(M*N);
score=zeros(1,numfids);
score2=zeros(1,numfids);
for z=1:numfids
    fid=fids(z);
    energies=zeros(1,runcheck);
    for g=1:runcheck
        var=load(['results/' filesave num2str(round(fid)) num2str(g)]);
        energies(g)=var.energy;
    end
    [~,I]=min(energies);
    var=load(['results/' filesave num2str(round(fid)) num2str(I)]);
    mapall=var.mapall;
    dict=var.dict;
    for i=1:total
          score(z)=score(z)+CI(i)*b2bmetric(truebetaEBSD(i,:),var.dict(var.mapall(i),:))^2;
          score2(z)=score2(z)+CI(i)*b2bmetric(dictp(mapallp(i),:),var.dict(var.mapall(i),:))^2;
    end
    score(z)=sqrt(score(z)/(M*N));
    score2(z)=sqrt(score2(z)/(M*N));
end
fid=fids(find(diff(score)>-5e-4,1));

save(['results/check' filesave],'score','score2','mapallp','dictp','kappap','truebetaEBSD')
% for z=1:numfids
%     fid=fids(z);
%     for g=1:runcheck
%         delete(['results/' filesave num2str(round(fid)) num2str(g) '.mat']);
%     end
% end

EBSDtemp=load(['../data/' filename 'EBSD.mat']);

EBSD=EBSDtemp.EBSD;
CI=EBSDtemp.CI;
[M,N]=size(CI);
beta=logical(EBSDtemp.betas);
codegenzaitzeff(M,N,0);
nr=ceil(M/40*gs/50);
nc=ceil(N/40*gs/50);
name=['results/' filesave num2str(round(fid))];

smallK=ceil((nr*nc)/8);
if numpar>1
    parpool(numpar)
    
    parfor pari=1:num
        tic;
        MStd(EBSD,CI,beta,fid,filesave,dt,dx,dy,nr,nc,pari);
        timings(pari)=toc;
    end
    [I,conval,conmap]=confidencemap(name,M,N,smallK,num,numpar);
    
    poolobj = gcp('nocreate');
    delete(poolobj);
    
else
    for i=1:num
        tic;
        MStd(EBSD,CI,beta,fid,filesave,dt,dx,dy,nr,nc,i);
        timings(i)=toc;
    end
    [I,conval,conmap]=confidencemap(name,M,N,smallK,num,numpar);
end

vars=load(['results/' filesave num2str(round(fid)) num2str(I)]);



mapall=vars.mapall;
dict=vars.dict;
energy=vars.energy;


betaEBSD=zeros(M*N,4);
K = size(dict,1);
alpha=reshape(EBSD, [M*N,3]);
alpha=E313toq(alpha);

[~,~,~,coords,sizecoords,~]=  bndcoords(mapall,K);
for k=1:K
    mu=dict(k,:);
    indices=coords(k,1:sizecoords(k));
    mask=beta(indices);
    acoord=indices(~mask);
    [best,~]=alpbbeta(alpha(acoord,:),mu);
    betaEBSD(acoord,:)=best;
    bcoord=indices(mask);
    betaEBSD(bcoord,:)=alpha(bcoord,:);

    
end
betaEBSD=qtoE313(betaEBSD);
betaEBSD=reshape(betaEBSD,[M,N,3]);
save(['results/' filesave num2str(round(fid))],'mapall','betaEBSD','dict','energy','conval','conmap','timings');

for i=1:num
    delete(['results/' filesave num2str(round(fid)) num2str(i) '.mat']);
end

