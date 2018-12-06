function MStd(EBSD,CI,beta,fid,Ks,filesave,dx,dy,dt,num,enec)

if nargin<7
    dx=1/100;
    dy=1/100;
end
if nargin<8
    dt=.02;
end
if nargin<9
    num=1;
end

if nargin<10
    enec=0;
end
[m,n,~]=size(EBSD);
dtstop=2^(-10);
try


        

    [mapall,dict,kappa]=initializeEBSDfast(EBSD,CI,beta,Ks,20,50,50);
    [mapall,dict,kappa]=regionmerging(mapall,dict,kappa,sEBSD,sCI,Ks,.5);

        

    %[mapall,dict,kappa,energy]=EBSDMStdfast(mapall,EBSD,CI,beta,dict,kappa,fid,dx*cs,dy*cs,dt/(4^z),ceil(20/cs),dtstop,mexed);
    [mapall,dict,kappa]=EBSDMStd(mapall,sEBSD,sCI,dict,kappa,fid,dt,dx,dy,dtstop,enec);

    energy=EBSDtdE(mapall,sEBSD,sCI,dict,fid,dtstop,dx,dy,enec);
    %energy=EBSDtdEfast(mapall,EBSD,CI,dict,fid,dx,dy,mexed);

catch
    mapall=zeros(m,n);
    dict={};
    energy=inf;
end
save(['results/' filesave num2str(fid) num2str(num)],'mapall','dict','energy');
end
