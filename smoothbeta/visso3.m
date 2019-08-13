function [im, oM]=visso3(u)
% visualizes an so3-valued image as an image of the corresponding euler
% angles
% [c,d,m,n] = size(u);
% assert(c == 3 & d == 3);
% 
% im = zeros(m,n,3);
% for i = 1:m
%     for j = 1:n
%         im(i,j,:) = rotm2eul(u(:,:,i,j));
%     end
% end
% imagesc(im(:,:,1)); colorbar

% imagesc(squeeze(u(1,1,:,:))); colorbar

addpath('../anglelib/')
%addpath('../../mtex-5.1.1/')

[c,d,m,n] = size(u);
assert(c == 3 & d == 3);

bungeimage = permute(reshape(BungeOfRMat(reshape(u, [3,3,m*n]), 'radians'), [3,m,n]), [2 3 1]);

[im, oM] = crystalcolormaps(bungeimage, 0);
%imagesc(im);

end