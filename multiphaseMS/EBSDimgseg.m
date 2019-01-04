function [mapall,dict,energy,gsizes]=EBSDimgseg(EBSD,CI,beta,fid,dt,dx,dy,nr,nc)
    dtstop=2^(-10);
    between=2;
    nt=6;
    [mapall,dict,kappa,mid]=initializeEBSDfast(EBSD,CI,beta,nr,nc);
    [mapall,dict,kappa]=regionmerging(mapall,dict,kappa,nr,nc,mid,fid);
    [mapall,dict,~,energy,gsizes]=EBSDMStdfast(mapall,EBSD,CI,beta,dict,kappa,fid,dt,dx,dy,dtstop,nt,between);
end