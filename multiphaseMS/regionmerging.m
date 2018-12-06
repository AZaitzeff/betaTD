function [mapall,newdict,newkappa]=regionmerging(mapall,dict,kappa,Ks,thres)


if nargin<7
    thres=3;
end
thres=thres*pi/360;%degrees to radians.
T=alphatobetatrans();
Pm=getsymmetries('cubic');
Pall=zeros(4,4,144);
for i=1:144
    Pall(:,:,i)=T(:,:,mod(i-1,6)+1)'*Pm(:,:,ceil(i/6));
end
%[m,n]=size(mapall);%size of level set function
% [~,~,z]=size(EBSD);
% EBSDflat=reshape(EBSD,[m*n,z]);
% EBSDflat=E313toq(EBSDflat);
% CIflat=reshape(CI,[m*n,1]);

m=Ks(2);
n=Ks(1);
K=m*n;
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
            values(z)=b2bmetric(dict(ind1,:),dict(ind2,:));
            ind(1,z)=ind1;
            ind(2,z)=ind2;
            z=z+1;
        end
        if i<m
            ind2=sub2ind([m,n],i+1,j);
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
        val=b2bmetric(dict(r1,:),dict(r2,:));
        if val<thres
            map(r2)=r1;
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



