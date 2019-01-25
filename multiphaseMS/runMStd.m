function runMStd(filename,filesave,numpar,num,runcheck,thres)

timings=zeros(1,num);
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
dt=2^-5;
nr=16;
nc=16;
%[mapallp,dictp,kappap,~]=initializeEBSDfast_mex(EBSD,CI,beta,nr,nc);
%truebetaEBSD=converttobetamap(EBSD,beta,dictp,mapallp);

next=[[7,8,9,10,11];[14,16,18,20,22];[29,33,37,41,45];[58,66,74,82,90];[118,136,154,172,190];...
    [236,272,308,344,380];[500,600,700,800,900]];

dts=[2^-3,2^-3,2^-4,2^-4,2^-5,2^-5];
fids=[12,25,50,100,200,400];
for iter=1:2
    
numfids=numel(fids);

totalcheck=runcheck*numfids;
if numpar>1
    parpool([1 numpar])
    
    parfor pari=1:totalcheck
        fidz=fids(ceil(pari/runcheck));
        dt=dts(ceil(pari/runcheck));
        MStd(EBSD,CI,beta,fidz,filesave,dt,dx,dy,nr,nc,mod(pari-1,runcheck)+1);
    end
    
    poolobj = gcp('nocreate');
    delete(poolobj);
    
else
    for i=1:totalcheck
        fid=fids(ceil(i/runcheck));
        dt=dts(ceil(i/runcheck));
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
    score(z)=prctile(vals,1);
    %mapall=var.mapall;
    gsizes(z)=round(prctile(var.gsizes,5));
    %save(['results/' filesave 'iter' num2str(round(fid))],'mapall');
end


if iter==1
I=find(score<thres,1)-1;
if isempty(I)
    I=numfids;
    startfid=fids(I);
    gsstart=gsizes(I);
    dt=2^-6;
elseif I==0
    startfid=6;
    dt=1/4;
    gsstart=20;
else
    startfid=fids(I);
    dt=dts(I);
    gsstart=gsizes(I);
end
dts(:)=dt;
fids=next(I+1,:);

else
    I=find(score<thres,1)-1;
    if isempty(I)
        fid=fids(numfids);
        I=numfids;
        gs=gsizes(I);
    elseif I==0
        fid=startfid;
        gs=gsstart;

    else
        fid=fids(I);
        gs=gsizes(I);
    end
end

end




for z=1:numfids
    tempfid=fids(z);
    for g=1:runcheck
        delete(['results/' filesave num2str(round(tempfid)) num2str(g) '.mat']);
    end
end


for tempfid=[12,25,50,100,200,400]
    for g=1:runcheck
        delete(['results/' filesave num2str(round(tempfid)) num2str(g) '.mat']);
    end
end



EBSDtemp=load(['../data/' filename 'EBSD.mat']);

EBSD=EBSDtemp.EBSD;
CI=EBSDtemp.CI;
[M,N]=size(CI);
beta=logical(EBSDtemp.betas);
codegenzaitzeff(M,N);
nr=ceil(M/sqrt(gs));
nc=ceil(N/sqrt(gs));
name=['results/' filesave num2str(round(fid))];

if numpar>1
    parpool([1 numpar])
    
    parfor pari=1:num
        tic;
        MStd(EBSD,CI,beta,fid,filesave,dt,dx,dy,nr,nc,pari);
        timings(pari)=toc;
    end
    
    poolobj = gcp('nocreate');
    delete(poolobj);
    
else
    for i=1:num
        tic;
        MStd(EBSD,CI,beta,fid,filesave,dt,dx,dy,nr,nc,i);
        timings(i)=toc;
    end
    
end
[I,conval,conmap]=confidencemap2(name,M,N,num);
[~,bndval,bndmap]=probmetric(name,w,num);
vars=load(['results/' filesave num2str(round(fid)) num2str(I)]);





mapall=vars.mapall;
dict=vars.dict;
energy=vars.energy;

betaEBSD=converttobetamap(EBSD,beta,dict,mapall);

save(['results/' filesave 'beta'],'mapall','betaEBSD','dict','energy','conval','conmap','bndval','bndmap','fid','score','gs','timings');

for i=1:num
    delete(['results/' filesave num2str(round(fid)) num2str(i) '.mat']);
end

