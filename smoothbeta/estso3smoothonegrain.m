function [newu,boxindex,fullsize]=estso3smoothonegrain(M,N,indices,EBSDflat,beta,CI,preindices,smoothEBSD,nt,dt,lam)                                

[ms,ns]=ind2sub([M,N],indices);
minysm=min(ms);
maxysm=max(ms);
minxsm=min(ns);
maxxsm=max(ns);
msm=maxysm-minysm+1;
nsm=maxxsm-minxsm+1;
uin=zeros(3,3,msm,nsm);
wsm=zeros(msm,nsm);

smlinind=sub2ind([msm,nsm],ms-minysm+1,ns-minxsm+1);
flatquat=zeros(4,msm*nsm);
betasm=zeros(1,msm*nsm);
backupindex=uint32(1);

[cols,rows]=meshgrid(1:nsm,1:msm);
boxindex=sub2ind([M,N],rows(:)+minysm-1,cols(:)+minxsm-1);
indexmap=buildindexmap(boxindex,preindices);
for i=1:msm
    for j=1:nsm
        index=sub2ind([M,N],i+minysm-1,j+minxsm-1);
        smindex=sub2ind([msm,nsm],i,j);
        flatquat(:,smindex)=EBSDflat(:,index);
        betasm(smindex)=beta(index);
        tempindex=indexmap(smindex);    
        uin(:,:,i,j)= RMatOfQuat(smoothEBSD(:,tempindex));

    end
end
wsm(smlinind)=CI(indices);

Rmat=bestparentmatrix(msm,nsm,flatquat,betasm,uin,0);
[u,~,~]=so3implicitfid(nt,dt,Rmat,uin,lam, wsm);
fullsize=msm*nsm;
newu=zeros(4,fullsize);
for smz=1:fullsize
    [ism,jsm]=ind2sub([msm,nsm],smz);
    newu(:,smz)=QuatOfRMat(u(:,:,ism,jsm));
end
                                