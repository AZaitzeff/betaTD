function [mapall,dict,energy,gsizes,flag]=EBSDimgseg(EBSD,CI,beta,fid,dt,dx,dy,nr,nc,numsub)
    dtstop=2^(-10);
    between=2;
    nt=6;
    [mapall,dict,kappa,mid]=initializeEBSDfast(EBSD,CI,beta,nr,nc,numsub);
    [mapall,dict,kappa]=regionmerging(mapall,dict,kappa,nr,nc,mid,fid);
    [mapall,dict,~,energy,gsizes,flag]=EBSDMStdfast(mapall,EBSD,CI,beta,dict,kappa,fid,dt,dx,dy,dtstop,nt,between,numsub);

end