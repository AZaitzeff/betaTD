function MStd(EBSD,CI,beta,fid,Ks,filesave,dt,dx,dy,num)

if nargin<7
    dx=1/100;
    dy=1/100;
end
if nargin<8
    dt=2^(-6);
end
if nargin<9
    num=1;
end

%try

    [mapall,dict,energy]=EBSDimgseg_mex(EBSD,CI,beta,fid,Ks,dt,dx,dy);

%catch
%    mapall=zeros(m,n);
%    dict={};
%    energy=inf;
%end
save(['results/' filesave num2str(round(fid)) num2str(num)],'mapall','dict','energy');
end
