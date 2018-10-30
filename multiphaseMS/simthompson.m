function x=simthompson(dt,dx,z,bound)
b=2/dt+2/dx^2;
a=-1/dx^2;
n=numel(z);
x=zeros(n,1);
if n==1
   x(1)=z(1)/b;
   return
end
c=zeros(n,1);
d=zeros(n,1);

if bound==1 || bound==3
    c(1)=2*a/(b); 
else
    c(1)=a/b;
end
d(1)=z(1)/(b);

for i=2:n-1
    c(i)=a/(b-a*c(i-1));
    d(i)=(z(i)-a*d(i-1))/(b-a*c(i-1));
end
if bound==2 || bound==3
    x(n)=(z(n)-2*a*d(n-1))/(b-2*a*c(n-1));
else
    x(n)=(z(n)-a*d(n-1))/(b-a*c(n-1));
end
for i=n-1:-1:1
    x(i)=d(i)-c(i)*x(i+1);
end
end