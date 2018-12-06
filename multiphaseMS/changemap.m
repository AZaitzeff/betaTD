function mapall=changemap(mapall,map)
[M,N]=size(mapall);
for i=1:M
    for j=1:N
        mapall(i,j)=map(mapall(i,j));      
    end
end

