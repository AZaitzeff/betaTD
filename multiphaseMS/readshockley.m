function E=readshockley(b1,b2,enec)

theta=b2bmetric(b1,b2)*2;

if enec>0
    b=enec^2;

    c=(1-b)*6/pi;
    E=sqrt(theta*c+b);

else
    b=exp(abs(enec));
    c=(exp(1)-b)*6/pi;
    E=log(theta*c+b);
end
% if theta<(1/180*pi)
%  E=2*pi*1/180*pi*(1-log(2*pi*1/180*pi));
% elseif theta<1/(2*pi)
%  E=2*pi*theta*(1-log(2*pi*theta));
% else
%  E=1;
% end
