function E=readshockley(b1,b2)

theta=b2bmetric(b1,b2)*2;
 if theta<(1/180*pi)
  E=1/4;
 else 
  E=(theta)*90/pi*1/2;
 end
% if theta<(1/180*pi)
%  E=2*pi*1/180*pi*(1-log(2*pi*1/180*pi));
% elseif theta<1/(2*pi)
%  E=2*pi*theta*(1-log(2*pi*theta));
% else
%  E=1;
% end
