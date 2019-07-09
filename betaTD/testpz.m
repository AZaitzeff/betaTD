addpath('../anglelib/')
load(['../data/' 'AFone' 'EBSD'])
 [M,N]=size(IQ);
 codegenzaitzeffpilchank(M,N)
 CI=CIfunc(CI);
 beta=logical(betas);

 [mapall,dict,newmap,newbetadict,energy,gsizes,flag]=EBSDimgsegpz_mex(EBSD,IQ,CI,beta,75,2^(-4),1/100,1/100,200,3,5,3,3,1);