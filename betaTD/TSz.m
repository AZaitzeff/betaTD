function u=TSz(u,T,nt,dx,dy,xdir,ydir,xsizes,ysizes,m,n,dir)
dx2=dx*dx;
dy2=dy*dy;
dt=T/nt;
tempu=u;

N=size(xsizes,1);
for iters=1:nt
    
for xkey=1:N
    xn=xsizes(xkey);
    if xn==0
        break
    end
    i=xdir(xkey,1);
    xcoor=xdir(xkey,2:xn+1);
    z=u(i,xcoor)/dt;
    
    if dir==0
        hor=(xcoor(1)==1)+(xcoor(xn)==n)*2;
    else
        hor=3;
    end
    if hor==1
        z(xn)=z(xn)+u(i,xcoor(xn)+1)/dx2;
    elseif hor==2
        z(1)=z(1)+u(i,xcoor(1)-1)/dx2;
    elseif hor==0
        z(1)=z(1)+u(i,xcoor(1)-1)/dx2;
        z(xn)=z(xn)+u(i,xcoor(xn)+1)/dx2; 
    end


    tempu(i,xcoor)=simthompson(2*dt,dx,z,hor);
end
u=tempu;
M=size(ysizes,1);
for ykey=1:M
    ym=ysizes(ykey);
    if ym==0
        break
    end
    j=ydir(ykey,1);
    ycoor=ydir(ykey,2:ym+1);
    z=u(ycoor,j)/dt;
    
    if dir==0
        ver=(ycoor(1)==1)+(ycoor(ym)==m)*2;
    else
        ver=3;
    end
    
    if ver==1
        z(ym)=z(ym)+u(ycoor(ym)+1,j)/dy2;
    elseif ver==2
        z(1)=z(1)+u(ycoor(1)-1,j)/dy2;
    elseif ver==0
        z(1)=z(1)+u(ycoor(1)-1,j)/dy2;
        z(ym)=z(ym)+u(ycoor(ym)+1,j)/dy2;
    end

    tempu(ycoor,j)=simthompson(2*dt,dy,z,ver);
end
u=tempu;

end