function [mapall,dict,kappa,mid]=initializeEBSDfast(EBSD,CI,beta,nr,nc,numsub)

%subn=Ks(1)*10;
%subm=Ks(2)*10;
[m,n,z]=size(EBSD);
rd=ceil(m/nr);
cd=ceil(n/nc);
K=nr*nc;
mid=(rd+cd)/2;
mapall=zeros(m,n);
rows=linspace(0,m,nr+1);
rspace=rows(2)-rows(1);
rows(2:end-1)=rows(2:end-1)+(rand(1,nr-1)-.5)*rspace/2;
rows=round(rows);
cols=linspace(0,n,nc+1);
cspace=cols(2)-cols(1);
cols(2:end-1)=cols(2:end-1)+(rand(1,nc-1)-.5)*cspace/2;
cols=round(cols);
for i=1:nr
    row=(rows(i)+1):rows(i+1);
    for j=1:nc
        col=(cols(j)+1):cols(j+1);
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

[dict,kappa]=estimatebetasfast(EBSDflat,CIflat,beta,mapall,K,numsub);