function [I,conval,conmap]=confidencemapmin(name,M,N,K,num,numpar)
addpath('../anglelib/')
Pm=getsymmetries('cubic');
conmap=ones(M,N);
mapall=zeros(M*N,num);
curnum=num*K;
dictall=zeros(num*K,4);
start=0;
enevec=zeros(1,num);
for z=1:num
    vars=load([name num2str(z)]);
    enevec(z)=vars.energy;
    dict=vars.dict;
    spots=size(dict,1);
    if (start+spots)>curnum
        curnum=ceil(curnum*num/z+spots);
        dictall=growbetas(dictall,curnum,start);
    end
    dictall((start+1):(start+spots),:)=dict;
    mapall(:,z)=vars.mapall(:)+start;
    start=start+spots;
end
[~,I]=min(enevec);
total=M*N;
if numpar>1
    parfor ind=1:total
        minbeta=dictall(mapall(ind,I),:);
        betas=dictall(mapall(ind,:),:);
        conmap(ind)=sum(b2bmetric(betas,minbeta).^2)/(num-1);
        %conmap(ind)=(kap);
    end
else
    for ind=1:total
        minbeta=dictall(mapall(ind,I),:);
        betas=dictall(mapall(ind,:),:);
        conmap(ind)=sum(b2bmetric(betas,minbeta).^2)/(num-1);
        %conmap(ind)=(kap);
    end
end
conval=sqrt(total/sum(conmap(:)));
end
function [newdictall]=growbetas(dictall,curnum,start)
    newdictall=zeros(curnum,4);
    for ind=1:start
        newdictall(ind,:)=dictall(ind,:);
    end
end