function map=MStd(EBSD,CI,fid,Ks,filename,section,DT)


if nargin<8
    DT=.02;
end
[m,n,~]=size(EBSD);
%maxener=inf;
%bu={};
%bK=0;
%bdict=0;
num=size(Ks,2);
enevec=zeros(1,num);
maps=cell(1,num);
dicts=cell(1,num);
for i=1:num
    K=Ks{i};
%for ztest =1:max_check
[u,newdict,newkappa]=initializeEBSD(EBSD,CI,K);
[u,dict,~]=EBSDMStd(u,EBSD,CI,newdict,newkappa,fid,DT);
energy=EBSDtdE(u,EBSD,CI,dict,fid);
enevec(i)=energy;
%if energy<maxener
%    maxener=energy;
%    bu=u;
%    bdict=dict;
%end
%end
total=size(u,2);
map=zeros(m,n);
for z=1:total
    map=map+(u{z}>0)*z;
end
maps{i}=map;
dicts{i}=dict;
end
[~,I]=min(enevec);

map=maps{I};
dict=dicts{I};
save(['results/' filename 'part' num2str(section) num2str(fid) '.mat'],'map','dict');