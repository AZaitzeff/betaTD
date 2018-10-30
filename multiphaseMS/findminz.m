function [r,I]=findminz(a,b)
n=numel(a);
r=zeros(1,n);
I=zeros(1,n,'logical');
for i=1:n
    if a(i)<=b(i)
        r(i)=a(i);
        I(i)=0;
    else
        r(i)=b(i);
        I(i)=1;
    end
end
