function MStd(EBSD,CI,beta,fid,filesave,dt,dx,dy,nr,nc,num,part)


[mapall,dict,energy,gsizes]=EBSDimgseg_mex(EBSD,CI,beta,fid,dt,dx,dy,nr,nc,part);
%[mapall,dict,energy,gsizes]=EBSDimgseg(EBSD,CI,beta,fid,dt,dx,dy,nr,nc);
%on sim avg time C vs MAT

save(['results/' filesave num2str(round(fid)) num2str(num)],'mapall','dict','energy','gsizes','-v7.3');
end
