function u=ADIz(u,T,dx,dy,xdir,ydir,xsizes,ysizes,m,n)
dx2=dx*dx;
dy2=dy*dy;
nt=ceil(log2(T/dx^2));
dt=T/2^nt;
tempu=u;
for iters=1:2*nt
    
for xkey=1:4000
    xn=xsizes(xkey);
    if xn==0
        break
    end
    i=xdir(xkey,1);
    xcoor=xdir(xkey,2:xn+1);
    if i==1
        z=u(i,xcoor)*(-2/dx2+2/dt)+(u(i+1,xcoor)+u(i+1,xcoor))/dx2;
    elseif i==m
        z=u(i,xcoor)*(-2/dx2+2/dt)+(u(i-1,xcoor)+u(i-1,xcoor))/dx2;
    else
        z=u(i,xcoor)*(-2/dx2+2/dt)+(u(i+1,xcoor)+u(i-1,xcoor))/dx2;
    end
    hor=(xcoor(1)==1)+(xcoor(xn)==n)*2;
    if hor==1
        z(xn)=z(xn)+u(i,xcoor(xn)+1)/dx2;
    elseif hor==2
        z(1)=z(1)+u(i,xcoor(1)-1)/dx2;
    elseif hor==0
        z(1)=z(1)+u(i,xcoor(1)-1)/dx2;
        z(xn)=z(xn)+u(i,xcoor(xn)+1)/dx2; 
    end


    tempu(i,xcoor)=simthompson(dt,dx,z,hor);
end
u=tempu;

for ykey=1:4000
    ym=ysizes(ykey);
    if ym==0
        break
    end
    j=ydir(ykey,1);
    ycoor=ydir(ykey,2:ym+1);
    if j==1
        z=u(ycoor,j)*(-2/dy2+2/dt)+(u(ycoor,j+1)+u(ycoor,j+1))/dy2;
    elseif j==n
        z=u(ycoor,j)*(-2/dy2+2/dt)+(u(ycoor,j-1)+u(ycoor,j-1))/dy2;
    else
        z=u(ycoor,j)*(-2/dy2+2/dt)+(u(ycoor,j-1)+u(ycoor,j+1))/dy2;
    end
    ver=(ycoor(1)==1)+(ycoor(ym)==m)*2;
    if ver==1
        z(ym)=z(ym)+u(ycoor(ym)+1,j)/dy2;
    elseif ver==2
        z(1)=z(1)+u(ycoor(1)-1,j)/dy2;
    elseif ver==0
        z(1)=z(1)+u(ycoor(1)-1,j)/dy2;
        z(ym)=z(ym)+u(ycoor(ym)+1,j)/dy2;
    end

    tempu(ycoor,j)=simthompson(dt,dy,z,ver);
end
u=tempu;

if iters>2 && mod(iters,2)==0
    dt=dt*2;
end
end