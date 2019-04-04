function [mapall,newdict,newkappa]=regionmergingmem(mapall,dict,kappa,m,n,mid,fid)
%fid=150;
%constant=0.28211;
gs=1/100;
threshold=prctile(kappa,50)/2;
K=m*n;
curK=K;
active=K;
contains=zeros(K,20);
neighbors=zeros(K,20);
contwidth=[K,20];
neiwidth=[K,20];
map=1:K;
point=map;
contains(1:K,1)=map;
contsize=ones(1,K);
neisize=zeros(1,K);
total=2*(n-1)*(m-1)+m-1+n-1;
values=zeros(1,total);
ind=zeros(2,total);
z=1;
for j=1:n
    for i=1:m
        ind1=sub2ind([m,n],i,j);
        if i<m
            ind2=sub2ind([m,n],i+1,j);
            neisize(ind1)=neisize(ind1)+1;
            neighbors(ind1,neisize(ind1))=ind2;
            neisize(ind2)=neisize(ind2)+1;
            neighbors(ind2,neisize(ind2))=ind1;
            values(z)=b2bmetric(dict(ind1,:),dict(ind2,:));
            ind(1,z)=ind1;
            ind(2,z)=ind2;
            z=z+1;
        end
        if j<n
            ind2=sub2ind([m,n],i,j+1);
            neisize(ind1)=neisize(ind1)+1;
            neighbors(ind1,neisize(ind1))=ind2;
            neisize(ind2)=neisize(ind2)+1;
            neighbors(ind2,neisize(ind2))=ind1;
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
    x1=findroot(map,ind1);
    x2=findroot(map,ind2);
    
    if x1~=x2
        r1=point(x1);
        r2=point(x2);
        if kappa(x1)>threshold && kappa(x2)>threshold
            %per=numel(intersect(contains(r1,1:contsize(r1)),neighbors(r2,1:neisize(r2))));
            per=andsum(contains(r1,1:contsize(r1)),neighbors(r2,1:neisize(r2)),contsize(r1),neisize(r2));
           
            area1=contsize(r1);
            area2=contsize(r2);
            area=min(area1,area2);
            perterm=per*gs*mid;
            val=b2bmetric(dict(x1,:),dict(x2,:));
            fidterm=val*fid*area*(gs*mid)^2;
            if fidterm<=perterm
                if area1>=area2
                    combinearray(x1,x2);
                else
                    combinearray(x2,x1);
                end
                active=active-1;
            end
        else
            if kappa(x1)>=kappa(x2)
                combinearray(x1,x2);
            else
                combinearray(x2,x1);
            end
            active=active-1;
        end
    end
    if curK*.9>active
        simplify();
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

    function combinearray(x1,x2)
        map(x2)=x1;
        l1=point(x1);
        l2=point(x2);
        C = union(contains(l1,1:contsize(l1)),contains(l2,1:contsize(l2)));
        sizeC=numel(C);
        if sizeC>contwidth(2)
            contwidth(2)=ceil(1.1*sizeC);
            [contains,contsize]=growarray(contains,contsize,curK,contwidth);
        end
        contains(l1,1:sizeC)=C;
        contsize(l1)=sizeC;
        
        C = union(neighbors(l1,1:neisize(l1)),neighbors(l2,1:neisize(l2)));
        sizeC=numel(C);
        if sizeC>neiwidth(2)
            neiwidth(2)=ceil(1.1*sizeC);
            [neighbors,neisize]=growarray(neighbors,neisize,curK,neiwidth);
        end
        neighbors(l1,1:sizeC)=C;
        neisize(l1)=sizeC;
        
        
    end

    function simplify()
        newK=0;
        for kt=1:K
            if map(kt)==kt
                
                cadk=point(kt);
                newK=newK+1;
                point(kt)=newK;
                tempsize=contsize(cadk);
                contains(newK,1:tempsize)=contains(cadk,1:tempsize);
                contsize(newK)=tempsize;
                
                tempsize=neisize(cadk);
                neighbors(newK,1:tempsize)=neighbors(cadk,1:tempsize);
                neisize(newK)=tempsize;
            end
        end
        contains((newK+1):curK,:)=[];
        contsize((newK+1):curK)=[];
        neighbors((newK+1):curK,:)=[];
        neisize((newK+1):curK)=[];
        
        curK=newK;
        contwidth(1)=curK;
        neiwidth(1)=curK;
    end
end


function thesum=andsum(A,B,sizeA,sizeB)
thesum=0;
i=1;
j=1;
while i<=sizeA && j<=sizeB
    if A(i)<B(j)
        i=i+1;
    elseif A(i)>B(j)
        j=j+1;
    else
        thesum=thesum+1;
        i=i+1;
        j=j+1;
    end
end


end