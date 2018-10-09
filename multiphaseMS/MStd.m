function map=MStd(EBSD,CI,fid,fixeds,filename,section,max_check,DT)


if nargin<7
    max_check=1;
end
if nargin<8
    DT=.02;
end
[m,n,~]=size(EBSD);
maxener=inf;
bu={};
bK=0;
bdict=0;
num=size(fixeds,2);
for i=1:num
    fixed=fixeds(i,:);
for ztest =1:max_check
K=prod(fixed);
[u,newdict,newkappa]=initializeEBSD(EBSD,CI,K,fixed);
[u,dict,~]=EBSDMStd(u,EBSD,CI,newdict,newkappa,fid,DT);
energy=EBSDtdE(u,EBSD,CI,dict,fid);
if energy<maxener
    maxener=energy;
    bu=u;
    bdict=dict;
    bK=K;
end

end
end
map=zeros(m,n);
for i=1:bK
    map=map+(bu{i}>0)*i;
end
dict=bdict;
save(['results/' filename 'part' num2str(section) num2str(fid) '.mat'],'map','dict');