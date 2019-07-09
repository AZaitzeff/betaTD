function E=readshockley(b1,b2,enec)

theta=b2bmetric(b1,b2)*2;

if enec>0 && enec<.75
    b=enec^2;

    c=(1-b)*6/pi;
    E=sqrt(theta*c+b);

elseif enec<=0
    b=exp(abs(enec));
    c=(exp(1)-b)*6/pi;
    E=log(theta*c+b);
else
    c=enec;
     if theta<(1/360*pi)
        E=2*pi*1/360*pi/c*(1-log(2*pi*1/360*pi/c));
     elseif theta<1/(2*pi/c)
        E=2*pi*theta/c*(1-log(2*pi*theta/c));
     else
        E=1;
     end
end
% if theta<(1/180*pi)
%  E=2*pi*1/180*pi*(1-log(2*pi*1/180*pi));
% elseif theta<1/(2*pi)
%  E=2*pi*theta*(1-log(2*pi*theta));
% else
%  E=1;
% end
