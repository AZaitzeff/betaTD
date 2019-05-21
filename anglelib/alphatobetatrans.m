function T=alphatobetatrans()
    invD = [0.5406 0.7046 -0.0600 0.4558]'/norm([0.5406 0.7046 0.0600 0.4558]);
    %hexsym= [ 0, 0, 0, 1.;
    %v               0.,0.,-1/2,sqrt(3)/2;
    %                0.,0.,-sqrt(3)/2,1/2;
    %                1., 0., 0., 0;
    %                1/2,sqrt(3)/2,0.,0.;
    %                sqrt(3)/2,1/2,0.,0.;];
    hexsym=HexSymmetries();
    T=zeros(4,4,6);
    ImpIndexes = [2 3 9 11 8 6];
    %ImpIndexes = [1 7 9 8 5 6];
    for i=1:6
        T(:,:,i)=quatmatrix2(multquat(hexsym(:,ImpIndexes(i)),invD));
    end
end