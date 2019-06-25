function mapall=changemapleavezeros(mapall,map)
[M,N]=size(mapall);
for i=1:M
    for j=1:N
        val=mapall(i,j);
        if val~=0
            mapall(i,j)=map(val);
        end
    end
end