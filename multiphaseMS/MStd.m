function MStd(EBSD,CI,fid,Ks,filesave,dx,dy,dt,step,num,enec)

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
    enec=0;
end
[m,n,~]=size(EBSD);
try

    nt=log2(step);
    for z=0:nt
        cs=step/2^z;
         
        sEBSD=EBSD(1:cs:end,1:cs:end,:);
        sCI=CI(1:cs:end,1:cs:end,:);
        
        if z==0
            [mapall,dict,kappa]=initializeEBSDfast(sEBSD,sCI,Ks,20,50,50);
            [mapall,dict,kappa]=regionmerging(mapall,dict,kappa,sEBSD,sCI,Ks,.5);
        else
            [m,n]=size(sCI);
            mapall = imresize(mapall, [m n], 'nearest');
        end
        
        if z==nt
            dtstop=2^(-10);
        else
            dtstop=dt/(4^(z+2));
        end
        %[mapall,dict,kappa]=EBSDMStdfast(mapall,sEBSD,sCI,dict,kappa,fid,dx*cs,dy*cs,dt/(4^z),ceil(20/cs),dtstop,mexed);
        [mapall,dict,kappa]=EBSDMStd(mapall,sEBSD,sCI,dict,kappa,fid,dt,dx,dy,dtstop,enec);
    end
    energy=EBSDtdE(mapall,sEBSD,sCI,dict,fid,dtstop,dx,dy,enec);
    %energy=EBSDtdEfast(mapall,EBSD,CI,dict,fid,dx,dy,mexed);

catch
    mapall=zeros(m,n);
    dict={};
    energy=inf;
end
save(['results/' filesave num2str(fid) num2str(num)],'mapall','dict','energy');
end
