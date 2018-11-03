function [EBSD,CI]=cleanEBSDdata(EBSD,CI,CIthres,winlength)
if nargin<4
    winlength=2;
end
mask=CI<CIthres;
[M,N]=size(CI);
[row,col]=find(mask);
n=sum(mask(:));
for i=1:n
    y=max(row(i)-winlength,1):min(row(i)+winlength,M);
    x=max(col(i)-winlength,1):min(col(i)+winlength,N);
    tempCI=CI(y,x);
    [~,I]=max(tempCI(:));
    [r,c]=ind2sub([numel(y) numel(x)], I);
    temp=EBSD(y,x,:);
    EBSD(row(i),col(i),:)=temp(r,c,:);   
end
CI(mask)=0;