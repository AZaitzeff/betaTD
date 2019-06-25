function vals=datasamplez(ind,num,weights)
vals=zeros(num,1);
prob=cumsum(weights);
Z=numel(weights);
prob=prob/prob(Z);
for i=1:num
    p=rand();
    v=binserach(p,prob,floor((Z+1)/2),1,Z);
    
   vals(i)=ind(v); 
end
end
function ind=binserach(p,prob,z,l,u)
    if p>=prob(z)
        if p<prob(z+1)
            ind=z+1;
        else
            ind=binserach(p,prob,floor((z+u)/2),z,u);
        end
    else
        if z==l
            ind=z;
        else
            ind=binserach(p,prob,floor((z+l)/2),l,z);
        end
    end

end

