function [vals2,whichn]=matchmetric(mapall,dict)
%[M,N,z]=size(EBSD);
%EBSDflat=reshape(EBSD,[M*N,z]);
%EBSDflat=E313toq(EBSDflat);
[K,~]=size(dict);
%[~,~,~,coords,sizecoords,~]=  bndcoords(mapall,K);
[neighbors]=  findneigh(mapall,K);
%vals1=zeros(1,K);
vals2=zeros(1,0);
whichn=zeros(2,0);
z=1;
for i =1:K
    for j=1:K
        if neighbors(i,j)
%             val=0;
%             linind=coords(j,1:sizecoords(j));
%             bmask=beta(linind);
%             alphacoord=linind(~bmask);
%             betacoord=linind(bmask);
%             if ~isempty(alphacoord)
%                 val=val+sum(CI(alphacoord).*alpbmetric(EBSDflat(alphacoord,:),dict(i,:)));
%             end
%             if ~isempty(betacoord)
%                 val=val+sum(CI(betacoord).*b2bmetric(EBSDflat(betacoord,:),dict(i,:)));
%             end
%             whichn(1,z)=i;
%             whichn(2,z)=j;
%             vals1(z)=val/sum(CI(linind));
            vals2(z)=b2bmetric(dict(j,:),dict(i,:));
            z=z+1;
        end
    end
end

% for k=1:K
%     val=0;
%     linind=coords(k,1:sizecoords(k));
%     bmask=beta(linind);
%     alphacoord=linind(~bmask);
%     betacoord=linind(bmask);
%     if ~isempty(alphacoord)
%         val=val+sum(CI(alphacoord).*alpbmetric(EBSDflat(alphacoord,:),dict(k,:)));
%     end
%     if ~isempty(betacoord)
%         val=val+sum(CI(betacoord).*b2bmetric(EBSDflat(betacoord,:),dict(k,:)));
%     end
%     vals1(k)=val/sum(CI(linind));
%     
% end
%vals1=vals1*2*180/pi;
vals2=vals2*2*180/pi;