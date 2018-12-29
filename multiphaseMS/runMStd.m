function runMStd(filename,filesave,numpar,num,checknoise,gs)
if nargin<5
    checknoise=1;
end

if nargin<6
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
%codegenzaitzeff(M,N,1);
alpha=reshape(EBSD, [M*N,3]);
alpha=E313toq(alpha);
%betas=EBSDtemp.betas(rows,cols);
dt=2^-5;

%nr=ceil(M/10*gs/50);
%nc=ceil(N/10*gs/50);
%[mapallp,dictp,kappap,~]=initializeEBSDfast_mex(EBSD,CI,beta,nr,nc);
%truebetaEBSD=converttobetamap(EBSD,beta,dictp,mapallp);

nr=ceil(M/40*gs/50);
nc=ceil(N/40*gs/50);

fids=25:25:300;
numfids=numel(fids);
runcheck=3;
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
for z=1:numfids
    fid=fids(z);
    energies=zeros(1,runcheck);
    for g=1:runcheck
        var=load(['results/' filesave num2str(round(fid)) num2str(g)]);
        energies(g)=var.energy;
    end
    [~,I]=min(energies);
    var=load(['results/' filesave num2str(round(fid)) num2str(I)]);
    for i=1:total
          score(z)=score(z)+CI(i)*alpbmetric(alpha(i,:),var.dict(var.mapall(i),:))^2;
    end
    score(z)=sqrt(score(z)/(M*N));
end

fid=kneedle(fids,score);

save(['results/check' filesave],'score')
if checknoise
name=['results/' filesave num2str(round(fid))];
smallK=ceil((nr*nc)/8);
if numpar>1
    parpool(numpar)
    
    parfor pari=1:num
        MStd(EBSD,CI,beta,fid,filesave,dt,dx,dy,nr,nc,pari);
    end
    [I,conval,conmap]=confidencemap(name,M,N,smallK,num,numpar);
    
    poolobj = gcp('nocreate');
    delete(poolobj);
    
else
    for i=1:num
        MStd(EBSD,CI,beta,fid,filesave,dt,dx,dy,nr,nc,i);
    end
    [I,conval,conmap]=confidencemap(name,M,N,smallK,num,numpar);
end
end

if ~checknoise || conval<.05
    message='Ran on full dataset';

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
        MStd(EBSD,CI,beta,fid,filesave,dt,dx,dy,nr,nc,pari);
    end
    [I,conval,conmap]=confidencemap(name,M,N,smallK,num,numpar);
    
    poolobj = gcp('nocreate');
    delete(poolobj);
    
else
    for i=1:num
        MStd(EBSD,CI,beta,fid,filesave,dt,dx,dy,nr,nc,i);
    end
    [I,conval,conmap]=confidencemap(name,M,N,smallK,num,numpar);
end

vars=load(['results/' filesave num2str(round(fid)) num2str(I)]);



mapall=vars.mapall;
dict=vars.dict;
energy=vars.energy;

else
    message='dataset is too noisy, if you want to run full run again with checknoise=0';
end
betaEBSD=converttobetamap(EBSD,beta,dict,mapall);

save(['results/' filesave num2str(round(fid))],'mapall','betaEBSD','dict','energy','conval','conmap','message');

for i=1:num
    delete(['results/' filesave num2str(round(fid)) num2str(i) '.mat']);
end

