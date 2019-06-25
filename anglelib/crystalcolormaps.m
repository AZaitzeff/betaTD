function [colors,oM]=crystalcolormaps(EBSD,alpha)
    if alpha
         cs=loadCIF('Ti-Titanium-alpha.cif');
    else
        cs=loadCIF('Ti-Titanium-beta.cif');
    end
    if ndims(EBSD)==3
        [m,n,~]=size(EBSD);

        oM = ipfHSVKey(cs);
        ori = orientation('Euler',EBSD(:,:,1),EBSD(:,:,2),EBSD(:,:,3));
        colors = oM.orientation2color(ori);
        colors=reshape(colors,[m,n,3]);
    else
        oM = ipfHSVKey(cs);
        ori = orientation('Euler',EBSD(:,1),EBSD(:,2),EBSD(:,3));
        colors = oM.orientation2color(ori);

    end

end