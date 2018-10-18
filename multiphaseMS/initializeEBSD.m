function [u,newdict,newkappa]=initializeEBSD(EBSD,CI,Ks)

h=1/100;
[m,n,~]=size(EBSD);
x=(1:n)*h;
y=(1:m)*h;
[X,Y]=meshgrid(x,y);
xs=zeros(m*n,2);
xs(:,1)=X(:);
xs(:,2)=Y(:);
K=prod(Ks);
if numel(Ks)==2
    center=zeros(K,2);
    [xpts,ypts]=meshgrid(((1:Ks(1))-1/2)/Ks(1)*n*h,((1:Ks(2))-1/2)/Ks(2)*m*h);
    center(:,1)=xpts(:)+randn(K,1)*n*h/(Ks(1)*8);
    center(:,2)=ypts(:)+randn(K,1)*m*h/(Ks(2)*8);
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
newdict=cell(1,K);
newkappa=cell(1,K);
u=cell(1,K);
%arr=[1,2,4,5,6];
for i=1:K
    temp=val;
    temp(:,i)=[];
    u{i}=(reshape(min(temp,[],2)-val(:,i), [m n]));
    %newdict{i}=dict(arr(i));
    newdict{i}=dict(i);
    newkappa{i}=kappa(i);
end