function [u,newdict,newkappa]=initializeEBSD(EBSD,CI,K,fixed)
if nargin<4
    fixed=[0,0];
end
h=1/100;
[m,n,z]=size(EBSD);
x=(1:n)*h;
y=(1:m)*h;
[X,Y]=meshgrid(x,y);
xs=zeros(m*n,2);
xs(:,1)=X(:);
xs(:,2)=Y(:);
if prod(fixed)==K
    center=zeros(K,2);
    [xpts,ypts]=meshgrid(((1:fixed(1))-1/2)/fixed(1)*n*h,((1:fixed(2))-1/2)/fixed(2)*m*h);
    center(:,1)=xpts(:);
    center(:,2)=ypts(:);
else
Db=0;
for i=1:10
xc=rand(K,2);
xc(:,1)=xc(:,1)*n*h;
xc(:,2)=xc(:,2)*m*h;
D=pdist(xc);
Dm=squareform(D)+max(D(:))*eye(K);
Dtotal = sum(min(Dm));
if Dtotal>Db
    Db=Dtotal;
    center=xc;
end
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
[dict,kappa]=estimatebetas(EBSD,CI,zeros(m,n),mapall,[],0,1,200,16);
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