function map=MStd(EBSD,CI,fid,Ks,filename,section,max_check)


if nargin<7
    max_check=1;
end

[m,n,~]=size(EBSD);
maxener=inf;
bu={};
bK=0;
bdict=0;
for K=Ks
for ztest =1:max_check

[u,newdict,newkappa]=initializeEBSD(EBSD,CI,K);
[u,dict,~]=EBSDMStd(u,EBSD,CI,newdict,newkappa,fid);
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