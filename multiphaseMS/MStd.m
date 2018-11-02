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
try

    nt=log2(step);
    for z=0:nt
        cs=step/2^z;
         
        sEBSD=EBSD(1:cs:end,1:cs:end,:);
        sCI=CI(1:cs:end,1:cs:end,:);
        
        if z==0
            [mapall,dict,kappa]=initializeEBSDfast(sEBSD,sCI,Ks,16,50,50,mexed);
        else
            [m,n]=size(sCI);
            mapall = imresize(mapall, [m n], 'nearest');
        end
        
        if z==nt
            dtstop=2^(-12);
        else
            dtstop=dt/(2^z)/8;
        end
        [mapall,dict,kappa]=EBSDMStdfast(mapall,sEBSD,sCI,dict,kappa,fid,dx*cs,dy*cs,dt/(2^z),ceil(20/cs),dtstop,mexed);
    end
    energy=EBSDtdEfast(mapall,EBSD,CI,dict,fid,dx,dy,mexed);

catch
    mapall=zeros(m,n);
    dict={};
    energy=inf;
end
save(['results/' filesave num2str(fid) num2str(num)],'mapall','dict','energy');
end
