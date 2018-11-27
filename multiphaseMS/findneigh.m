function [neighbors]=  findneigh(mapall,K)
neighbors=zeros(K,K);
[M,N]=size(mapall);

for j=2:N-1
    for i =2:M-1
        k=mapall(i,j);
        neighbors(k,mapall(i+1,j))=1;
        neighbors(k,mapall(i,j+1))=1;
    end
end
neighbors=(neighbors+neighbors')>.5;
for k=1:K
    neighbors(k,k)=0;
end
end