function [mapall,dict,kappa,mid]=initializeEBSDfast(EBSD,CI,beta,nr,nc)

%subn=Ks(1)*10;
%subm=Ks(2)*10;
[m,n,z]=size(EBSD);
rd=ceil(m/nr);
cd=ceil(n/nc);
K=nr*nc;
mid=(rd+cd)/2;
mapall=zeros(m,n);
for i=1:nr
    for j=1:nc
        row=(i-1)*rd+1:min(i*rd,m);
        col=(j-1)*cd+1:min(j*cd,n);
        mapall(row,col)=i+(j-1)*nr;
        
    end
end


%x=(1:subn);
%y=(1:subm);
%[X,Y]=meshgrid(x,y);
%total=subn*subm;
%xs=zeros(subn*subm,2);
%xs(:,1)=X(:);
%xs(:,2)=Y(:);
%K=prod(Ks);
%center=zeros(K,2);
%[xpts,ypts]=meshgrid(((1:Ks(1))-1/2)/Ks(1)*subn,((1:Ks(2))-1/2)/Ks(2)*subm);
%center(:,1)=xpts(:)+randn(K,1)*subn/(Ks(1)*4);
%center(:,2)=ypts(:)+randn(K,1)*subm/(Ks(2)*4);

%val=zeros(subm*subn,K);
%for k=1:K
%    for i=1:total
%        val(i,k)=sum((xs(i,:)-center(k,:)).^2);
%    end
%end
%[~,r]=min(val,[],2);
%smallmap=reshape(r,[subm subn]);


%[m,n,~]=size(EBSD);
%x=1:Ks(1);
%y=1:Ks(2);
%K=prod(Ks);
%smallmap=(x-1)*Ks(2)+y';
%mapall = imresize(smallmap, [m n], 'nearest');
EBSDflat=reshape(EBSD, [m*n,z]);
EBSDflat=E313toq(EBSDflat);
CIflat=reshape(CI, [m*n,1]);

[dict,kappa]=estimatebetasfast(EBSDflat,CIflat,beta,mapall,K);