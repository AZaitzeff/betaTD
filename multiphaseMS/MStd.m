function MStd(EBSD,CI,fid,Ks,filesave,dx,dy,dt,step,num,reestbeta)

if nargin<7
    dx=1/100;
    dy=1/100;
end
if nargin<8
    dt=.02;
end
if nargin<9
    step=1;
end
if nargin<10
    num=1;
end
if nargin<11
    reestbeta=0;
end
[m,n,~]=size(EBSD);

if step>1
    sEBSD=EBSD(1:step:end,1:step:end,:);
    sCI=CI(1:step:end,1:step:end,:);
end

enevec=zeros(1,num);
maps=cell(1,num);
dicts=cell(1,num);
parfor i=1:num
if step==0
    [mapall,newdict,newkappa]=initializeEBSD(EBSD,CI,Ks);
    [mapall,dict,~]=EBSDMStd(mapall,EBSD,CI,newdict,newkappa,fid,dx,dy,dt);
    energy=EBSDtdE(mapall,EBSD,CI,dict,fid,dx,dy);
    enevec(i)=energy;
elseif step<0
    [mapall,newdict,newkappa]=initializeEBSDfast(EBSD,CI,Ks);
    [mapall,dict,~]=EBSDMStd(mapall,EBSD,CI,newdict,newkappa,fid,dx,dy,dt,10);
    energy=EBSDtdEfast(mapall,EBSD,CI,dict,fid,dx,dy);
    enevec(i)=energy;
    
else
    [smallmap,newdict,newkappa]=initializeEBSD(sEBSD,sCI,Ks);
    [smallmap,dict,kappa]=EBSDMStd(smallmap,sEBSD,sCI,newdict,newkappa,fid,dx*step,dy*step,dt,(2^-12)*step);
    mapall = imresize(smallmap, [m n], 'nearest');
    dt2=dt/step;
    [mapall,dict,~]=EBSDMStdfast(mapall,EBSD,CI,dict,kappa,fid,dx,dy,dt2);
    energy=EBSDtdEfast(mapall,EBSD,CI,dict,fid,dx,dy);
    enevec(i)=energy;
end
maps{i}=mapall;
dicts{i}=dict;
end
[~,I]=min(enevec);

mapall=maps{I};
dict=dicts{I};
if reestbeta
    [dict]=estimatebetas(EBSD,CI,zeros(size(CI)),mapall,[],0,1,500);
end
betaEBSD=zeros(m*n,4);
[clusterlist,~,~] = unique(mapall);
alpha=reshape(EBSD, [m*n,z]);
alpha=E313toq(alpha);
for z=clusterlist'
    mu=dict(z);
    zfun=find(z==mapall);

    [best,~]=alpbbeta(alpha(zfun,:),mu);
    betaEBSD(zfun,:)=best;

end
betaEBSD=qtoE313(betaEBSD);
betaEBSD=reshape(betaEBSD,[m,n,3]);
save(['results/' filesave num2str(fid)],'mapall','betaEBSD','dict');