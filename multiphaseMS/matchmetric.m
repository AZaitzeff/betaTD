function [vals2,whichn]=matchmetric(mapall,dict)

[K,~]=size(dict);
[neighbors]=  findneigh(mapall,K);
vals2=zeros(1,0);
whichn=zeros(2,0);
z=1;
for i =1:K
    for j=1:K
        if neighbors(i,j)
            vals2(z)=b2bmetric(dict(j,:),dict(i,:));
            z=z+1;
        end
    end
end


vals2=vals2*2*180/pi;