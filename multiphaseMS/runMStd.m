function runMStd(filename,filesave,numpar,num,runcheck,checknoise)
if nargin<5
    checknoise=1;
end

addpath('../anglelib/')


if nargin<3
    numpar=10;
    num=10;
end


EBSDtemp=load(['../data/' filename 'EBSD.mat']);
addpath('../anglelib/')
dx=1/100;
dy=dx*EBSDtemp.scale;

[M,N]=size(EBSDtemp.CI);
midm=ceil(M/2);
midn=ceil(N/2);
lm=max(midm-200,1);
um=min(midm+200,M);
ln=max(midn-200,1);
un=min(midn+200,N);
suby=lm:um;
subx=ln:un;
M=um-lm+1;
N=un-ln+1;
EBSD=EBSDtemp.EBSD(suby,subx,:);
CI=EBSDtemp.CI(suby,subx);
beta=logical(EBSDtemp.betas(suby,subx));
codegenzaitzeff(M,N);
%alpha=reshape(EBSD, [M*N,3]);
%alpha=E313toq(alpha);
%betas=EBSDtemp.betas(rows,cols);
%dts=[2^-5 2^-5.33 2^-5.66 2^-6];
dt=2^-4;
nr=20;
nc=20;
%[mapallp,dictp,kappap,~]=initializeEBSDfast_mex(EBSD,CI,beta,nr,nc);
%truebetaEBSD=converttobetamap(EBSD,beta,dictp,mapallp);

totalfid=6;
allfids=zeros(1,totalfid);
allscore=zeros(1,totalfid);
allgrains=zeros(1,totalfid);
len=100;
fids=[150,300];
numfids=numel(fids);
con=1;
while len>20
    for fid=fids
        allfids(con)=fid;
        con=con+1;
    end
    con=con-2;
totalcheck=runcheck*numfids;
if numpar>1
    parpool([1 numpar])
    
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


score=zeros(1,numfids);
gsizes=zeros(1,numfids);
for z=1:numfids
    fid=fids(z);
    energies=zeros(1,runcheck);
    for g=1:runcheck
        var=load(['results/' filesave num2str(round(fid)) num2str(g)]);
        energies(g)=var.energy;
    end
    [~,I]=min(energies);
    var=load(['results/' filesave num2str(round(fid)) num2str(I)]);
    [vals,~]=matchmetric(var.mapall,var.dict);
    score(z)=prctile(vals,2);
    %mapall=var.mapall;
    gsizes(z)=round(prctile(var.gsizes,5));
    %save(['results/' filesave 'iter' num2str(round(fid))],'mapall');
end

for iter=1:2
    allscore(con)=score(iter);
    allgrains(con)=gsizes(iter);
    con=con+1;
end

len=len/2;
if score(1)>=1.2 && score(2)>=1.2
    fids=[len+fids(2),2*len+fids(2)];
elseif score(1)<=.8 && score(2)<=.8
    fids=[-len+fids(1),-2*len+fids(1)];
elseif score(1)>=.8 && score(2)<=.8
    fids=[-len+fids(1),len+fids(1)];
elseif score(1)>=1.2 && score(2)<=1.2    
    fids=[-len+fids(2),len+fids(2)];
else
    fids=[len+fids(1),-len+fids(2)];
end
end

[fids,I]=sort(allfids);
tempscore=allscore(I);
score=zeros(1,totalfid);
for i=1:totalfid
    if i==1
        score(1)=3/4*tempscore(1)+tempscore(2)/4;
    elseif i==totalfid
        score(totalfid)=3/4*tempscore(totalfid)+tempscore(totalfid-1)/4;
        
    else
        score(i)=tempscore(i-1)/4+1/2*tempscore(i)+tempscore(i+1)/4;
    end
end
gsizes=allgrains(I);

I=find(score<1,1)-1;
if isempty(I)
    I=totalfid;
elseif I==0
    I=1;
end
fid=fids(I);
gs=gsizes(I);
if checknoise
    name=['results/' filesave num2str(round(fid))];
    smallK=ceil((nr*nc)/8);
    if numpar>1
        parpool([1 numpar])
        [I,conval,conmap]=confidencemap(name,M,N,smallK,num,numpar);
        poolobj = gcp('nocreate');
        delete(poolobj);

    else
        [I,conval,conmap]=confidencemap(name,M,N,smallK,num,numpar);
    end
    vars=load(['results/' filesave num2str(round(fid)) num2str(I)]);
end

for z=1:numfids
    tempfid=fids(z);
    for g=1:runcheck
        delete(['results/' filesave num2str(round(tempfid)) num2str(g) '.mat']);
    end
end


if ~checknoise || conval<.05
    message='Ran on full dataset';

EBSDtemp=load(['../data/' filename 'EBSD.mat']);

EBSD=EBSDtemp.EBSD;
CI=EBSDtemp.CI;
[M,N]=size(CI);
beta=logical(EBSDtemp.betas);
codegenzaitzeff(M,N);
nr=ceil(M/sqrt(gs));
nc=ceil(N/sqrt(gs));
name=['results/' filesave num2str(round(fid))];

smallK=ceil((nr*nc)/8);
if numpar>1
    parpool([1 numpar])
    
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

else
    message='dataset is too noisy, if you want to run on full, run again with checknoise=0';
end




mapall=vars.mapall;
dict=vars.dict;
energy=vars.energy;

betaEBSD=converttobetamap(EBSD,beta,dict,mapall);

save(['results/' filesave 'beta'],'mapall','betaEBSD','dict','energy','conval','conmap','fid','fids','score','gs','message');

for i=1:num
    delete(['results/' filesave num2str(round(fid)) num2str(i) '.mat']);
end

