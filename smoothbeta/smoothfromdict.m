function [smoothEBSD,K] = smoothfromdict(mapall,dict)
%SMOOTHFROMDICT Summary of this function goes here
%   Detailed explanation goes here
[z,K]=size(dict);

[M,N]=size(mapall);
smoothEBSD=zeros(z,M*N);

[~,~,~,coords,sizecoords,~]=  bndcoords(mapall,K);

for k=1:K
    smoothEBSD(1,coords(k,1:sizecoords(k)))=dict(1,k);
    smoothEBSD(2,coords(k,1:sizecoords(k)))=dict(2,k);
    smoothEBSD(3,coords(k,1:sizecoords(k)))=dict(3,k);
    smoothEBSD(4,coords(k,1:sizecoords(k)))=dict(4,k);
end

