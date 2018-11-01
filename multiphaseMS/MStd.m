function MStd(EBSD,CI,fid,Ks,filesave,dx,dy,dt,step,num,mexed)

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
    mexed=0;
end
[m,n,~]=size(EBSD);

if step==1
    [mapall,newdict,newkappa]=initializeEBSD(EBSD,CI,Ks);
    [mapall,dict,~]=EBSDMStd(mapall,EBSD,CI,newdict,newkappa,fid,dx,dy,dt);
    energy=EBSDtdE(mapall,EBSD,CI,dict,fid,dx,dy);
elseif step<1
    [mapall,newdict,newkappa]=initializeEBSDfast(EBSD,CI,Ks,16,50,50,mexed);
    [mapall,dict,~]=EBSDMStd(mapall,EBSD,CI,newdict,newkappa,fid,dx,dy,dt,20,mexed);
    energy=EBSDtdEfast(mapall,EBSD,CI,dict,fid,dx,dy,mexed);
    
else
    sEBSD=EBSD(1:step:end,1:step:end,:);
    sCI=CI(1:step:end,1:step:end,:);
    [smallmap,newdict,newkappa]=initializeEBSD(sEBSD,sCI,Ks);
    [smallmap,dict,kappa]=EBSDMStd(smallmap,sEBSD,sCI,newdict,newkappa,fid,dx*step,dy*step,dt,(2^-12)*step);
    mapall = imresize(smallmap, [m n], 'nearest');
    dt2=dt/(2*step);
    [mapall,dict,~]=EBSDMStdfast(mapall,EBSD,CI,dict,kappa,fid,dx,dy,dt2,0,mexed);
    energy=EBSDtdEfast(mapall,EBSD,CI,dict,fid,dx,dy,mexed);
end
save(['results/' filesave num2str(fid) num2str(num)],'mapall','dict','energy');
end
