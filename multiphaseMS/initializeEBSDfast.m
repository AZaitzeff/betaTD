function [mapall,dict,kappa]=initializeEBSDfast(EBSD,CI,Ks,check,subm,subn,mexed)
if nargin<4
    check=16;
end
if nargin<6
    subn=50;
    subm=50;
end
if nargin<7
    mexed=0;
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
[dict,kappa]=estimatebetasfast(EBSD,CI,mapall,check,mexed);