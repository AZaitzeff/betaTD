function [mapall,dict,kappa,K]=paircheck(mapall,dict,kappa,CIflat,EBSDflat,K,dt,dx,dy,enec)
T=alphatobetatrans();
Pm=getsymmetries('cubic');
Pall=zeros(4,4,144);
for i=1:144
    Pall(:,:,i)=T(:,:,mod(i-1,6)+1)'*Pm(:,:,ceil(i/6));
end

numsub=400;
thres=10;
[neighbors]=  findneigh(mapall,K);
map=1:K;
total=sum(neighbors(:))/2;
pairs=zeros(total,2);
vals=zeros(total,1);
z=1;
for k1=1:K
    for k2=(k1+1):K
        if neighbors(k1,k2)
            pairs(z,1)=k1;
            pairs(z,2)=k2;
            vals(z)=b2bmetric(dict(k1,:),dict(k2,:));
            z=z+1;
        end
    end
end
[newvals,I]=sort(vals);
pairs=pairs(I,:);
u=cell(K,1);
for k=1:K
    u{k}=1*(mapall==k);
end
[ls]=td2dz(u,dt,dx,dy);
for z=1:total
    if newvals>thres/360*pi
        break
    end
    k1=findroot(map,pairs(z,1));
    k2=findroot(map,pairs(z,2));
    if k1~=k2
        if enec
            energy=readshockley(dict(k1,:),dict(k2,:));
        else
            energy=1;
        end
        mask1=mapall==k1;
        mask2=mapall==k2;
        first=energy*(sum(ls{k1}(mask2))+sum(ls{k2}(mask1)))+...
            sum(CIflat(mask1(:)).*alpbmetric(EBSDflat(mask1(:),:),dict(k1,:))')+...
            sum(CIflat(mask2(:)).*alpbmetric(EBSDflat(mask2(:),:),dict(k2,:))');

        indices=find(mask1+mask2);
        newind=datasamplez(indices,numsub,CIflat(indices));
        EBSDtemp=EBSDflat(newind,:);
        [newg1, kap1, ~] = VMFEMfast(EBSDtemp, Pall,1,dict(k1,:),kappa(k1));
        [newg2, kap2, ~] = VMFEMfast(EBSDtemp, Pall,1,dict(k2,:),kappa(k2));
        second=sum(CIflat(indices).*alpbmetric(EBSDflat(indices,:),newg1)');
        third=sum(CIflat(indices).*alpbmetric(EBSDflat(indices,:),newg2)');
        if second<third
            val=second;
            newg=newg1;
            kap=kap1;
        else
            val=third;
            newg=newg1;
            kap=kap2;
        end
        if val<first || energy<(2*pi*1/180*pi*(1-log(2*pi*1/180*pi))+1e-4)
            dict(k1,:)=newg;
            kappa(k1)=kap;
            mapall(mapall==k2)=k1;
            map(k2)=k1;

        end
    end
inds=unique(map);
K=numel(inds);
newdict=zeros(K,4);
newkappa=zeros(K,1);
newmapall=mapall;

for k=1:K
    ind=inds(k);
    newkappa(k)=kappa(ind);
    newdict(k,:)=dict(ind,:);
    newmapall(mapall==ind)=k;
end
end
