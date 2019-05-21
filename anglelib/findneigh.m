function [neigharray,sizeneigh]=  findneigh(mapall,K)
neighbors=eye(K);
[M,N]=size(mapall);

for j=2:N-1
    for i =2:M-1
        k=mapall(i,j);
    	if k>0 && mapall(i+1,j)>0
        	neighbors(k,mapall(i+1,j))=1;
        end
        if k>0 && mapall(i,j+1)>0
        	neighbors(k,mapall(i,j+1))=1;
        end
    end
end
neighbors=(neighbors+neighbors')>.5;

sizeneigh=sum(neighbors)-1;
neigharray=zeros(K,max(sizeneigh));
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