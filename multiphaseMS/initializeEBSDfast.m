function [mapall,newdict,newkappa]=initializeEBSDfast(EBSD,CI,Ks,numsub,check,subm,subn)
if nargin<4
    numsub=200;
end
if nargin<5
    check=16;
end
if nargin<7
    subn=50;
    subm=50;
end
[m,n,~]=size(EBSD);

x=(1:subn);
y=(1:subm);
[X,Y]=meshgrid(x,y);
xs=zeros(subn*subm,2);
xs(:,1)=X(:);
xs(:,2)=Y(:);
K=prod(Ks);
center=zeros(K,2);
[xpts,ypts]=meshgrid(((1:Ks(1))-1/2)/Ks(1)*subn,((1:Ks(2))-1/2)/Ks(2)*subm);
center(:,1)=xpts(:)+randn(K,1)*subn/(Ks(1)*4);
center(:,2)=ypts(:)+randn(K,1)*subm/(Ks(2)*4);

val=zeros(subm*subn,K);
for k=1:K
    mean=sum((xs-center(k,:)).^2,2);
    val(:,k)=mean;
end
[~,r]=min(val,[],2);
smallmap=reshape(r,[subm subn]);
%[m,n,~]=size(EBSD);
%x=1:Ks(1);
%y=1:Ks(2);
%K=prod(Ks);
%smallmap=(x-1)*Ks(2)+y';
mapall = imresize(smallmap, [m n], 'nearest');
[dict,kappa]=estimatebetas(EBSD,CI,zeros(m,n),mapall,[],0,1,numsub,check);
newdict=cell(1,K);
newkappa=cell(1,K);

for i=1:K
    newdict{i}=dict(i);
    newkappa{i}=kappa(i);
end