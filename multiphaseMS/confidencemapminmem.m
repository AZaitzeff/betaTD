function [I,conval,conmap]=confidencemapminmem(name,M,N,num)
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
mainK=size(maindict,1);
[~,~,~,maincoords,mainsizecoords,~]=  bndcoords(vars.mapall,mainK);
total=M*N;
for z=1:num
    if z~=I
        vars=load([name num2str(z)]);
        dict=vars.dict;
        K=size(dict,1);
        [~,~,~,coords,sizecoords,~]=  bndcoords(vars.mapall,K);
        for i=1:mainK
            for j=1:K
                ind = intersect(maincoords(i,1:mainsizecoords(i)),coords(j,1:sizecoords(j)));
                conmap(ind)=conmap(ind)+b2bmetric(maindict(i,:),dict(j,:))^2;
            end
        end
            
    end
end
conmap=conmap/(num-1);
conval=sqrt(total/sum(conmap(:)));
end