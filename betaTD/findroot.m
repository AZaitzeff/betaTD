function ind=findroot(map,ind)
    tind=map(ind);
    while ind~=tind
        ind=tind;
        tind=map(ind);
    end
end