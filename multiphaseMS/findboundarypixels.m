function mask=findboundarypixels(u)
dy=diff(u(:,2:end-1),1,1);
dx=diff(u(2:end-1,:),1,2);

mask= (dy(1:end-1,:)==1) | (dy(2:end,:)==-1) | (dx(:,1:end-1)==1) | (dx(:,2:end)==-1);
mask=padarray(mask,[1,1]);