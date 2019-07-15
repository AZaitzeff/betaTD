function [neigharray,sizeneigh]=  findneigh(mapall,K)
neighbors=spalloc(K,K,K*25);
%neighbors=zeros(K,K);
[M,N]=size(mapall);

for j=2:N-1
    for i =2:M-1
        k=mapall(i,j);
        if k>0 && mapall(i+1,j)>0 && mapall(i+1,j)~=k
            neighbors(k,mapall(i+1,j))=1;
            neighbors(mapall(i+1,j),k)=1;
        end
        if k>0 && mapall(i,j+1)>0 && mapall(i,j+1)~=k
            neighbors(k,mapall(i,j+1))=1;
            neighbors(mapall(i,j+1),k)=1;
        end
    end
end

sizeneigh=full(sum(neighbors));
neigharray=zeros(K,max(sizeneigh(:)));
for k=1:K
    z=1;
    for k2=1:K
        if k2~=k && neighbors(k,k2)
            neigharray(k,z)=k2;
            z=z+1;
        end

    end
end
end