function [I,conval,conmap]=confidencemapminC(name,M,N,num)
addpath('../anglelib/')
conmap=zeros(M,N);
enevec=zeros(1,num);
maxK=0;
for z=1:num
    vars=load([name num2str(z)]);
    enevec(z)=vars.energy;
    maxK=max(maxK,size(vars.dict,1));
end
[~,I]=min(enevec);
vars=load([name num2str(I)]);
maindict=vars.dict;
mainmap=vars.mapall;
total=M*N;

vectorType1 = coder.typeof(1, [M N], [false false]);
vectorType3 = coder.typeof(1, [maxK 4], [true false]);

%vectorType3 = coder.typeof(1, [K 4], [true false]);
%vectorType4 = coder.typeof(1, [K 1], [true false]);
%codegen regionmerging -args {vectorType1,vectorType3,vectorType4,[1 1],1}
%codegen EBSDMStdfast -args {vectorType1,vectorType2,vectorType1,vectorType1l,vectorType3,vectorType4,1,1,1,1,1,1,1}
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
conval=sqrt(total/sum(conmap(:)));
end