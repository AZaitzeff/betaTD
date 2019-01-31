function runMStd(filename,filesave,numpar,num,runcheck)

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
codegenzaitzeff(M,N);
%alpha=reshape(EBSD, [M*N,3]);
%alpha=E313toq(alpha);
%betas=EBSDtemp.betas(rows,cols);
%dts=[2^-5 2^-5.33 2^-5.66 2^-6];
nr=ceil(M/20);
nc=ceil(N/20);
%[mapallp,dictp,kappap,~]=initializeEBSDfast_mex(EBSD,CI,beta,nr,nc);
%truebetaEBSD=converttobetamap(EBSD,beta,dictp,mapallp);

%next=[[7,8,9,10,11];[14,16,18,20,22];[29,33,37,41,45];[58,66,74,82,90];[118,136,154,172,190];...
%    [236,272,308,344,380];[500,550,600,650,700]];
%fac=[1/8,1/(4*sqrt(2)),1/4,1/(2*sqrt(2)),1/2,1/sqrt(2),1,sqrt(2),2,2*sqrt(2)];
dts=[2^-3,2^-3,2^-4,2^-4,2^-4,2^-4,2^-5,2^-5,2^-5,2^-5];
fids=[25,25*sqrt(2),50,50*sqrt(2),100,100*sqrt(2),200,200*sqrt(2),400,400*sqrt(2)];
    
numfids=numel(fids);

totalcheck=runcheck*numfids;
% if numpar>1
%     parpool([1 numpar])
%     
%     parfor pari=1:totalcheck
%         fidz=fids(ceil(pari/runcheck));
%         dt=dts(ceil(pari/runcheck));
%         MStd(EBSD,CI,beta,fidz,filesave,dt,dx,dy,nr,nc,mod(pari-1,runcheck)+1);
%     end
%     
%     poolobj = gcp('nocreate');
%     delete(poolobj);
%     
% else
%     for i=1:totalcheck
%         fid=fids(ceil(i/runcheck));
%         dt=dts(ceil(i/runcheck));
%         MStd(EBSD,CI,beta,fid,filesave,dt,dx,dy,nr,nc,mod(i-1,runcheck)+1);
%     end
% end

score=zeros(1,numfids);
flags=zeros(1,numfids,'logical');
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
    if var.flag
        [val,~]=fidmetric(EBSD,CI,beta,var.mapall,var.dict);
        score(z)=val;
        flags(z)=1;
    else
        flags(z)=0;
    end
    flags(z)=var.flag;
    %mapall=var.mapall;
    gsizes(z)=round(prctile(var.gsizes,5));
    %save(['results/' filesave 'iter' num2str(round(fid))],'mapall');
end

ind=bestbreaks(score(flags));
tempind=1:numfids;
tempind=tempind(flags);
ind=tempind(ind);
% confidence=zeros(1,2);
% %prob=zeros(1,2);
% if numpar>1
%     parpool([1 numpar])
% for i=1:2
%     fid=fids(ind(i));
%     name=['results/' filesave num2str(round(fid))];
%     [I,conval,~]=confidencemapmin(name,M,N,40,runcheck,numpar);
%     confidence(i)=conval;
% end
%     poolobj = gcp('nocreate');
%     delete(poolobj);
%     
% else
%     for i=1:2
%         fid=fids(ind(i));
%         name=['results/' filesave num2str(round(fid))];
%         [I,conval,~]=confidencemapmin(name,M,N,40,runcheck,numpar);
%         confidence(i)=conval;
%     end
% end
% 
% if confidence(2)>confidence(1)
%     fid=fids(ind(2));
%     gs=max(gsizes(ind(2))*1.2,400);
% else
%     fid=fids(ind(1));
%     gs=max(gsizes(ind(1))*1.2,400);
% end


fids(ind)
gsizes(ind)



for z=1:numfids
    tempfid=fids(z);
    for g=1:runcheck
        delete(['results/' filesave num2str(round(tempfid)) num2str(g) '.mat']);
    end
end


for iters=1:2
fid=fids(ind(iters));
gs=max(gsizes(ind(iters))*.9,324);
dt=dts(ind(iters));

EBSDtemp=load(['../data/' filename 'EBSD.mat']);

EBSD=EBSDtemp.EBSD;
CI=EBSDtemp.CI;
[M,N]=size(CI);
beta=logical(EBSDtemp.betas);
codegenzaitzeff(M,N);
nr=ceil(M/sqrt(gs));
nc=ceil(N/sqrt(gs));
smallK=ceil((nr*nc)/8);
name=['results/' filesave num2str(round(fid))];

if numpar>1
    parpool([1 numpar])
    
    parfor pari=1:num
        tic;
        MStd(EBSD,CI,beta,fid,filesave,dt,dx,dy,nr,nc,pari);
        timings(pari)=toc;
    end
    [I,conval,conmap]=confidencemapmin(name,M,N,smallK,num,numpar);
    poolobj = gcp('nocreate');
    delete(poolobj);
    
else
    for i=1:num
        tic;
        MStd(EBSD,CI,beta,fid,filesave,dt,dx,dy,nr,nc,i);
        timings(i)=toc;
    end
    [I,conval,conmap]=confidencemapmin(name,M,N,smallK,num,numpar);
    
end
w=5;
[~,bndval,bndmap]=probmetric(name,w,num);
vars=load(['results/' filesave num2str(round(fid)) num2str(I)]);



mapall=vars.mapall;
dict=vars.dict;
energy=vars.energy;

betaEBSD=converttobetamap(EBSD,beta,dict,mapall);

save(['results/' filesave 'beta' num2str(round(fid))],'mapall','betaEBSD','dict','energy','conval','conmap','bndval','bndmap','fid','score','gs','timings');

for i=1:num
    delete(['results/' filesave num2str(round(fid)) num2str(i) '.mat']);
end

end
