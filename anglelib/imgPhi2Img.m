function [Img] = imgPhi2Img(Phi, GoodVoxels, drt, isCubic)
    refDir = zeros(3,1);
    if(drt==1 || drt==2 || drt==3)
        refDir(round(drt))=1;
    else % set the reference direction as the mean direction
        Phi1=Phi(:,:,1);
        Phi2=Phi(:,:,2);
        Phi3=Phi(:,:,3);
        P = [Phi1(find(GoodVoxels)) Phi2(find(GoodVoxels)) Phi3(find(GoodVoxels))];
        refDir = mean(P);
    end
    R = EulerAng2IPF(Phi, refDir, isCubic);
    Img=zeros(size(R));
    Img(:,:,1) = R(:,:,1).*GoodVoxels;
    Img(:,:,2) = R(:,:,2).*GoodVoxels;
    Img(:,:,3) = R(:,:,3).*GoodVoxels;
end