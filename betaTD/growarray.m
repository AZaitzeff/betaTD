function [new,sizenew]=growarray(orig,sizeorig,os,ns,ind1,ind2)
if nargin<5
    ind1=0;
end
if nargin<6
    ind2=0;
end
new=zeros(ns);
if ~ind1
    sizenew=zeros(ns(1),1);
else
    sizenew=sizeorig;
end

for i =1:os
    v=sizeorig(i)+ind1;
    if ~ind1
        sizenew(i)=v;
    end
    if ind2
        new(i,1:v+1)=orig(i,1:v+1);
    else
        new(i,1:v)=orig(i,1:v);
    end
end
