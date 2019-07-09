function [map,dict,K]=cleanmap(map,dict,K,minsize)

[~,~,~,~,sizecoords,~]=  bndcoords(map,K);

replacemap=zeros(1,K);

z=1;
for k=1:K
    if sizecoords(k)>=minsize
        replacemap(k)=z;
        z=z+1;
        dict(:,z)=dict(:,k);
    end
end
map=changemapleavezeros(map,replacemap);
dict(:,z:K)=[];
K=z-1;
        
