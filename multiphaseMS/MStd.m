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
%try


        

    [mapall,dict,kappa]=initializeEBSDfast(EBSD,CI,beta,Ks,20,50,50,400);
    [mapall,dict,kappa]=regionmerging(mapall,dict,kappa,sEBSD,sCI,Ks,.5);

        

    %[mapall,dict,kappa,energy]=EBSDMStdfast_mex(mapall,EBSD,CI,beta,dict,kappa,fid,dx*cs,dy*cs,dt/(4^z),ceil(20/cs),dtstop,mexed);
    [mapall,dict,~]=EBSDMStd(mapall,EBSD,CI,dict,kappa,fid,dt,dx,dy,dtstop,enec);

    energy=EBSDtdE(mapall,EBSD,CI,dict,fid,dtstop,dx,dy,enec);
    %energy=EBSDtdEfast(mapall,EBSD,CI,dict,fid,dx,dy,mexed);

%catch
%    mapall=zeros(m,n);
%    dict={};
%    energy=inf;
%end
save(['results/' filesave num2str(fid) num2str(num)],'mapall','dict','energy');
end
