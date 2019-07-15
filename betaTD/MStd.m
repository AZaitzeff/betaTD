function MStd(EBSD,CI,beta,fid,filesave,dt,scale,nr,nc,num)
[mapall,dict,energy,gsizes]=EBSDimgseg_mex(EBSD,CI,beta,fid,dt,scale,nr,nc);
save(['results/' filesave num2str(round(fid)) num2str(num)],'mapall','dict','energy','gsizes','-v7.3');
end
