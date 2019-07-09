function [I,conval,conmap]=confidencemap2(name,M,N,num)
addpath('../anglelib/')

conmap=zeros(M,N);
total=M*N;
enevec=zeros(1,num);
for z=1:num
    vars=load([name num2str(z)]);
    enevec(z)=vars.energy;
end
[~,I]=min(enevec);
varmin=load([name num2str(I)]);
for z=1:num
    if z~=I
        vars=load([name num2str(z)]);
        for ind=1:total
            conmap(ind)=conmap(ind)+(b2bmetric(varmin.dict(varmin.mapall(ind),:),vars.dict(vars.mapall(ind),:)).^2);
        end
    end
end

conmap(ind)=conmap(ind)/(num-1);
        %conmap(ind)=(kap);

conval=sqrt(total/sum(conmap(:)));
