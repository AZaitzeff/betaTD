function LA=buildindexmap(A,B)
%this is naive make it better (current finds closest index instead of nearest row col)
n1=numel(A);
n2=numel(B);
LA=zeros(1,n1,'uint32');
z=1;
for i=1:n1
    while z<=n2 && B(z)<A(i)
        z=z+1;
    end
    if z>n2
        LA(i)=n2;
    elseif B(z)==A(i)
        LA(i)=z;
    elseif z>1 && abs(B(z)-A(i))>abs(B(z-1)-A(i))
        LA(i)=z-1;
    else
        LA(i)=z;
    end
end