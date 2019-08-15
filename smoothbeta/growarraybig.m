function [new]=growarraybig(orig,sizeorig,os,ns)

new=zeros(4,ns(2),ns(1));
for i =1:os
    v=sizeorig(i);
    new(:,1:v,i)=orig(:,1:v,i);
end