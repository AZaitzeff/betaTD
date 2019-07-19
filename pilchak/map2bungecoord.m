function [EBSD,valued]=map2bungecoord(map,dict,filled,K)
    [M,N]=size(map);
    EBSD=zeros(M,N,3);
    EBSD1=zeros(M,N);
    EBSD2=zeros(M,N);
    EBSD3=zeros(M,N);
    valued=zeros(M,N,'logical');
    [~,~,~,coords,sizecoords,~]=  bndcoords(map,K);
    for k =1:K
        if filled(k)
            EulerB=qtoE313(dict(:,k));
            kcoord=coords(k,1:sizecoords(k));
            EBSD1(kcoord)=EulerB(1);
            EBSD2(kcoord)=EulerB(2);
            EBSD3(kcoord)=EulerB(3);
            valued(kcoord)=1;
        end
    end
    EBSD(:,:,1)=EBSD1;
    EBSD(:,:,2)=EBSD2;
    EBSD(:,:,3)=EBSD3;
end