function conmap=calc2norms(conmap,maindict,mainmap,dict,map,total)
    for ind=1:total
        minbeta=maindict(:,mainmap(ind));
        betas=dict(:,map(ind));
        conmap(ind)=conmap(ind)+b2bmetric(betas,minbeta)^2;
    end
end