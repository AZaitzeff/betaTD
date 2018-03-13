function [gridu,g1,g2]=oneregion(region,EBSD,CI,betas,Img,bd,fid,uinb,g1,g2,plot)
    [n,m,~]=size(EBSD);
    area=sum(sum(region));
    if area>100
        [row,col]=find(region);
        yb=max(min(row)-7,1):min(max(row)+7,n);
        xb=max(min(col)-7,1):min(max(col)+7,m);
        z=region(yb,xb);
        sEBSD=EBSD(yb,xb,:);
        sCI=CI(yb,xb);
        sbetas=betas(yb,xb);
        se = strel('disk',5);
        dilatedreg=imdilate(z,se);
        uin=uinb(yb,xb);
        [u,g1,g2] = phiupdate(100000,1/(2*100^2),uin,z,dilatedreg,sEBSD,sCI,sbetas,fid,g1,g2);
        gridu=zeros(n,m);
        gridu(yb,xb)=u;
    else
        gridu=region;
    end
    if plot
        figure('visible','off');
        subplot(121)
        imshow(Img(yb,xb,:));
        subplot(122)
        imshow(1-bd(yb,xb));hold on;
        contour(u, [.5 .5], 'r');contour(z, [.5 .5], 'b');
        print(['test' num2str(fid)],'-dpng');
        close
    end