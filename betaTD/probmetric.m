function [I,score,probmap]=probmetric(name,w,num)
addpath('../anglelib/')
SE = strel('disk',w,4);
enevec=zeros(1,num);
for z=1:num
    vars=load([name num2str(z)]);
    enevec(z)=vars.energy;
end
[~,I]=min(enevec);
vars=load([name num2str(I)]);
[~,bnds]=  bndmap(vars.mapall,[],[]);
template = imdilate(bnds,SE);
probmap=template;
score=0;
for z=1:num
    if z~=I
        vars=load([name num2str(z)]);
        [~,bnds]=  bndmap(vars.mapall,[],[]);
        check = imdilate(bnds,SE);
        total=sum(check(:)|template(:));
        if total>0
            score=score+sum(check(:)&template(:))/ total;
        else
            score=score+1;
        end
        probmap=probmap+check;
    end
end

probmap=probmap/num;
score=score/(num-1);