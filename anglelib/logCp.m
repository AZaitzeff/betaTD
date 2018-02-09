function [y] = logCp(kappa, p)
%   This function is used for calculate the density of Von-Mises Fisher
%   distribution.  
% 
% inputs,
%   kappa : It is the kappa parameter of the Von-Mises Fisher distribution.
%
%   p : It is the dimension of Von-Mises Fisher distribution which is an
%       integer.
%
% outputs,
%   y : It is the output value of this function.
%
% Function is written by Yu-Hui Chen, University of Michigan

y = log((2*pi)^(-p/2) * (kappa.^(p/2-1) ./ besseli(p/2-1, kappa)));

KK = kappa;
KK(kappa<=700)=[];

y2=log((2*pi)^(-p/2) * (700.^(p/2-1) ./ besseli(p/2-1, 700)));
y1=log((2*pi)^(-p/2) * (699.^(p/2-1) ./ besseli(p/2-1, 699)));
yKK = y2+(y2-y1)*(KK-700);
y(kappa>700) = yKK;

end