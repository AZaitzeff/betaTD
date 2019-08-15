function plotsmoothwbnds(mapall,smoothEBSD)
% visualizes an q-valued image as an image of the corresponding euler
[M,N]=size(mapall);
[colorsbnd,bnds]=  bndmap(mapall,[],[]);
bungegrains=qtoE313(smoothEBSD);

bungegrains=reshape(bungegrains',[M,N,3]);
[im, ~] = crystalcolormaps(bungegrains, 0);

imagesc(im);hold on
imagesc(colorsbnd,'AlphaData',bnds)

%imagesc(im);

end