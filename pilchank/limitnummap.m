function [map,dict,K]=limitnummap(map,dict,K,numwanted)
if numwanted>K
    return
end
[~,~,~,~,sizecoords,~]=  bndcoords(map,K);
[~,ind]=sort(sizecoords,'descend');
replacemap=zeros(1,K);
newdict=zeros(4,numwanted);
for z=1:numwanted
    k=ind(z);
    replacemap(k)=z;
    newdict(:,z)=dict(:,k);
end
map=changemapleavezeros(map,replacemap);
dict=newdict;
K=numwanted;
end