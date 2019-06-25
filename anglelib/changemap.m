function mapall=changemap(mapall,map)
[M,N]=size(mapall);
for i=1:M
    for j=1:N
        val=map(mapall(i,j));
        if val~=0
            mapall(i,j)=val; 
        else
            if j<(N/2)
                ind=0;
                while val==0
                    ind=ind+1;
                    val=map(mapall(i,j+ind));
                end
            else
                val=mapall(i,j-1);
            end
            mapall(i,j)=val;  
            
		end 
    end
end

