function [I,conval,conmap]=confidencemapminC(name,M,N,num)
conmap=zeros(M,N);
enevec=zeros(1,num);
maxK=0;
for z=1:num
    vars=load([name num2str(z)]);
    enevec(z)=vars.energy;
    maxK=max(maxK,size(vars.dict,2));
end
[~,I]=min(enevec);
vars=load([name num2str(I)]);
maindict=vars.dict;
mainmap=vars.mapall;
total=M*N;

vectorType1 = coder.typeof(1, [M N], [false false]);
vectorType3 = coder.typeof(1, [4 maxK], [false true]);

codegen calc2norms -args {vectorType1,vectorType3,vectorType1,vectorType3,vectorType1,1}

for z=1:num
    if z~=I
        vars=load([name num2str(z)]);
        dict=vars.dict;
        map=vars.mapall;
        conmap=calc2norms_mex(conmap,maindict,mainmap,dict,map,total);
    end
end
conmap=conmap/(num-1);
conval=sqrt(sum(conmap(:))/total);
end