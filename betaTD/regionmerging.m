function [mapall,newdict,newkappa]=regionmerging(mapall,dict,kappa,m,n,mid,fid)
%fid=150;
%constant=0.28211;
gs=1/100;
threshold=prctile(kappa,50)/2;
K=m*n;
contains=eye(K,'logical');
neighbors=zeros(K,K,'logical');
map=1:K;
total=2*(n-1)*(m-1)+m-1+n-1;
values=zeros(1,total);
ind=zeros(2,total);
z=1;
for i=1:m
    for j=1:n
        ind1=sub2ind([m,n],i,j);
        if j<n
            ind2=sub2ind([m,n],i,j+1);
            neighbors(ind1,ind2)=1;
            neighbors(ind2,ind1)=1;
            values(z)=b2bmetric(dict(ind1,:),dict(ind2,:));
            ind(1,z)=ind1;
            ind(2,z)=ind2;
            z=z+1;
        end
        if i<m
            ind2=sub2ind([m,n],i+1,j);
            neighbors(ind1,ind2)=1;
            neighbors(ind2,ind1)=1;
            values(z)=b2bmetric(dict(ind1,:),dict(ind2,:));
            ind(1,z)=ind1;
            ind(2,z)=ind2;
            z=z+1;
        end
        
    end
end
[~,order]=sort(values);
ind=ind(:,order);

for z =1:total
    ind1=ind(1,z);
    ind2=ind(2,z);
    r1=findroot(map,ind1);
    r2=findroot(map,ind2);
    if r1~=r2
        if kappa(r1)>threshold && kappa(r2)>threshold
            per=andsum(contains(r1,:),neighbors(r2,:),K);
            area1=sum(contains(r1,:));
            area2=sum(contains(r2,:));
            area=min(area1,area2);
            perterm=per*gs*mid;
            val=b2bmetric(dict(r1,:),dict(r2,:));
            fidterm=val*fid*area*(gs*mid)^2;
            if fidterm<=perterm
                if area1>=area2
                    map(r2)=r1;
                    contains(r1,:)=or(contains(r1,:),contains(r2,:));
                    neighbors(r1,:)=or(neighbors(r1,:),neighbors(r2,:));
                else
                    map(r1)=r2;
                    contains(r2,:)=or(contains(r1,:),contains(r2,:));
                    neighbors(r2,:)=or(neighbors(r1,:),neighbors(r2,:));
                end

            end
        else
            if kappa(r1)>=kappa(r2)
                map(r2)=r1;
                contains(r1,:)=or(contains(r1,:),contains(r2,:));
                neighbors(r1,:)=or(neighbors(r1,:),neighbors(r2,:));
            else
                map(r1)=r2;
                contains(r2,:)=or(contains(r1,:),contains(r2,:));
                neighbors(r2,:)=or(neighbors(r1,:),neighbors(r2,:));
            end
            
        end
    end

end



current=zeros(1,K);
for k=1:K
    if map(k)==k
        current(k)=1;
    end
end
newK=sum(current);
newdict=zeros(newK,4);
newkappa=zeros(newK,1);
newk=1;
newmap=1:K;
for k=1:K
    if current(k)
        newmap(k)=newk;
        newdict(newk,:)=dict(k,:);
        newkappa(newk)=kappa(k);
        newk=newk+1;
    end
end

for k=1:K
    newmap(k)=newmap(findroot(map,k));
end

mapall=changemap(mapall,newmap);
end

function thesum=andsum(A,B,K)
thesum=0;
for k=1:K
    if A(k) && B(k)
        thesum=thesum+1;
    end
end


end

