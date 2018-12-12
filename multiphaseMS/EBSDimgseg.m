function [mapall,dict,energy]=EBSDimgseg(EBSD,CI,beta,fid,Ks,dt,dx,dy)
    subsample=400;
    dtstop=2^(-10);
    between=2;
    nt=6;
    [mapall,dict,kappa]=initializeEBSDfast(EBSD,CI,beta,Ks,20,subsample);
    [mapall,dict,kappa]=regionmerging(mapall,dict,kappa,Ks,.5);
    [mapall,dict,~,energy]=EBSDMStdfast(mapall,EBSD,CI,beta,dict,kappa,fid,dt,dx,dy,dtstop,nt,between,subsample);
end