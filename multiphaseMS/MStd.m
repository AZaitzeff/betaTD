function MStd(EBSD,CI,beta,fid,filesave,dt,dx,dy,nr,nc,num)

%try

    [mapall,dict,energy]=EBSDimgseg_mex(EBSD,CI,beta,fid,dt,dx,dy,nr,nc);

%catch
%    mapall=zeros(m,n);
%    dict={};
%    energy=inf;
%end
save(['results/' filesave num2str(round(fid)) num2str(num)],'mapall','dict','energy','-v7.3');
end
