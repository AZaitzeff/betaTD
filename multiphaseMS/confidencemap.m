function [I,conval,conmap]=confidencemap(name,M,N,K,num,numpar)
addpath('../anglelib/')
Pm=getsymmetries('cubic');
vectorType1 = coder.typeof(1, [num 4], [false false]);
vectorType2 = coder.typeof(1, [4 4 24], [false false false]);
codegen VMFEMfast -args {vectorType1,vectorType2,1,vectorType1,5}
conmap=zeros(M,N);
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
    parfor ind=total
        betas=dictall(mapall(ind,:),:);
        [meanbeta, kap, ~] = VMFEMfast_mex(betas, Pm,-num,betas,50);
        %conmap(ind)=sum(b2bmetric(betas,meanbeta).^2)/num;
        conmap(ind)=(1/kap);
    end
else
    for ind=total
        betas=dictall(mapall(ind,:),:);
        [meanbeta, kap, ~] = VMFEMfast_mex(betas, Pm,-num,betas,50);
        %conmap(ind)=sum(b2bmetric(betas,meanbeta).^2)/num;
        conmap(ind)=(1/kap);
    end
end
conval=sqrt(sum(conmap(:)))/(M*N);
end
function [newdictall]=growbetas(dictall,curnum,start)
    newdictall=zeros(curnum,4);
    for ind=1:start
        newdictall(ind,:)=dictall(ind,:);
    end
end


