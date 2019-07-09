filesave='AFone';
load(['../data/' filesave 'EBSD'])
[M,N]=size(CI);
midm=ceil(M/2);
midn=ceil(N/2);
lm=max(midm-50,1);
um=min(midm+50,M);
ln=max(midn-50,1);
un=min(midn+50,N);
suby=lm:um;
subx=ln:un;
M=um-lm+1;
N=un-ln+1;
EBSD=EBSD(suby,subx,:);
CI=CI(suby,subx);
beta=logical(betas(suby,subx));
fid=75;

[mapall,dict,energy,gsizes,flag]=EBSDimgseg(EBSD,CI,beta,fid,2^-4,1/100,1/100,6,6,400);
