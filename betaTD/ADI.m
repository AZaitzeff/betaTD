function u=ADI(u,T,dx,dy,xdir,ydir,m,n)
dx2=dx*dx;
dy2=dy*dy;
nt=ceil(log2(T/dx^2));
dt=T/2^nt;
tempu=u;
nx=size(xdir,2);
ny=size(ydir,2);
for iters=1:(nt+1)
for xkey=1:nx
    values=xdir{xkey};
    i=values{1};
    
    for k=1:values{2}
        xcoor=values{k+2};
        xn=size(xcoor,2);
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
        %tempu(i+1,xcoor+1)=trisolve(dt,dx,z);
    end
end
u=tempu;

for ykey=1:ny
    values=ydir{ykey};
    j=values{1};
    for k=1:values{2}
        ycoor=values{k+2};
        ym=size(ycoor,2);
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
end
u=tempu;
if iters>1
    dt=dt*2;
end
end