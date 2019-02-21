function [I,conval,conmap]=confidencemapminmem(name,M,N,num,numpar)
addpath('../anglelib/')
conmap=zeros(M,N);
enevec=zeros(1,num);
for z=1:num
    vars=load([name num2str(z)]);
    enevec(z)=vars.energy;
end
[~,I]=min(enevec);
vars=load([name num2str(I)]);
maindict=vars.dict;
mainmap=vars.mapall;
total=M*N;
for z=1:num
    if z~=I
        vars=load([name num2str(z)]);
        dict=vars.dict;
        map=vars.mapall;
        if numpar>1
            parfor ind=1:total
                minbeta=maindict(mainmap(ind),:);
                betas=dict(map(ind),:);
                conmap(ind)=conmap(ind)+b2bmetric(betas,minbeta)^2;
                %conmap(ind)=(kap);
            end
        else
            for ind=1:total
                minbeta=maindict(mainmap(ind),:);
                betas=dict(map(ind),:);
                conmap(ind)=conmap(ind)+b2bmetric(betas,minbeta)^2;
                %conmap(ind)=(kap);
            end
        end
    end
end
conmap=conmap/(num-1);
conval=sqrt(total/sum(conmap(:)));
end