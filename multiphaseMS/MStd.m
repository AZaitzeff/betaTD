function map=MStd(EBSD,CI,fid,Ks,filename,section,max_check)


if nargin<7
    max_check=1;
end

[m,n,z]=size(EBSD);
x=(1:n)/100;
y=(1:m)/100;
[X,Y]=meshgrid(x,y);
pos=zeros(m*n,2);
pos(:,1)=X(:);
pos(:,2)=Y(:);
EBSDflat=reshape(EBSD,[m*n,z]);
CIflat=reshape(CI,[m*n,1]);
EBSDflat=E313toq(EBSDflat);
maxener=inf;
bu={};
bK=0;
bdict=0;
kmeans_fid=.16;
for K=Ks
for ztest =1:max_check

[mu,~,CIdx]=kmeansEBSD(EBSDflat,CIflat,pos,kmeans_fid,K,3,5);
u={};
dict={};
for i=1:K
    u{i}=(reshape(CIdx==i, [m n])*2.)-1;
    dict{i}=mu(i,:);
end

[u,dict]=EBSDMStd(u,EBSD,CI,dict,fid);
energy=EBSDtdE(u,EBSD,CI,dict,fid);
if energy<maxener
    maxener=energy;
    bu=u;
    bK=K;
    bdict=dict;
end

end
end
map=zeros(m,n);
for i=1:bK
    map=map+(bu{i}>0)*i;
end
dict=bdict;
save(['results/' filename 'part' num2str(section) num2str(fid) '.mat'],'map','dict');