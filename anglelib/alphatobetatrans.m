function T=alphatobetatrans()
    phi1=3*pi/4;
    theta=pi/2;
    phi2=65*pi/36;
    D=zeros(4,1);
    D(1)=cos((phi1-phi2)/2)*sin(theta/2);
    D(2)=sin((phi1-phi2)/2)*sin(theta/2);
    D(3)=sin((phi1+phi2)/2)*cos(theta/2);
    D(4)=cos((phi1+phi2)/2)*cos(theta/2);
    hexsym= [ 0, 0, 0, 1.;
                    0.,0.,-1/2,sqrt(3)/2;
                    0.,0.,-sqrt(3)/2,1/2;
                    1., 0., 0., 0;
                    1/2,sqrt(3)/2,0.,0.;
                    sqrt(3)/2,1/2,0.,0.;];
    invD=invquat(D);
    T=zeros(4,4,6);
    for i=1:6
        T(:,:,i)=quatmatrix(multquat(invD,hexsym(i,:)'));
    end
end