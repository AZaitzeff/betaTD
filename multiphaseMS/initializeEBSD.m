function [u,newdict,newkappa]=initializeEBSD(EBSD,CI,K)
h=1/100;
[m,n,z]=size(EBSD);
x=(1:n)/100;
y=(1:m)/100;
[X,Y]=meshgrid(x,y);
xs=zeros(m*n,2);
xs(:,1)=X(:);
xs(:,2)=Y(:);
Db=0;
for i=1:10
xc=rand(K,2);
xc(:,1)=xc(:,1)*n/100;
xc(:,2)=xc(:,2)*m/100;
D=pdist(xc);
Dm=squareform(D)+max(D(:))*eye(K);
Dtotal = sum(min(Dm));
if Dtotal>Db
    Db=Dtotal;
    center=xc;
end
end
val=zeros(m*n,K);
for k=1:K
    mean=sum((xs-center(k,:)).^2,2);
    val(:,k)=mean;
end
[~,r]=min(val,[],2);
mapall=zeros(m,n);
for i=1:K
    mapall=mapall+(reshape(r,[m n])==i)*i;
end
[dict,kappa]=estimatebetas(EBSD,CI,zeros(m,n),mapall);
newdict={};
newkappa={};
%arr=[1,2,4,5,6];
for i=1:K
    temp=val;
    temp(:,i)=[];
    u{i}=(reshape(min(temp,[],2)-val(:,i), [m n]));
    %newdict{i}=dict(arr(i));
    newdict{i}=dict(i);
    newkappa{i}=kappa(i);
end