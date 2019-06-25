function T=alphatobetatrans()
    %invD = [0.5406 0.7046 -0.0600 0.4558]'/norm([0.5406 0.7046 0.0600 0.4558]);
    phi1=3*pi/4;
    theta=pi/2;
    phi2=65*pi/36;
    invD=zeros(4,1);
    invD(1)=cos((phi1+phi2)/2)*cos(theta/2);
    invD(2)=-cos((phi1-phi2)/2)*sin(theta/2);
    invD(3)=-sin((phi1-phi2)/2)*sin(theta/2);
    invD(4)=-sin((phi1+phi2)/2)*cos(theta/2);
    %hexsym= [ 0, 0, 0, 1.;
    %v               0.,0.,-1/2,sqrt(3)/2;
    %                0.,0.,-sqrt(3)/2,1/2;
    %                1., 0., 0., 0;
    %                1/2,sqrt(3)/2,0.,0.;
    %                sqrt(3)/2,1/2,0.,0.;];
    hexsym=HexSymmetries();
    T=zeros(4,4,6);
    %ImpIndexes = [2 3 9 11 8 6];
    ImpIndexes = [2 6 8 3 9 11];
    for i=1:6
        T(:,:,i)=quatmatrix2(multquat(hexsym(:,ImpIndexes(i)),invD));
    end
end