function [y] = Ap(x, p)
%   This function is used for calculate the density of Von-Mises Fisher
%   distribution. It is basically the ratio between two modified Bessel 
%   functions of the first kind. 
% 
% inputs,
%   x : It is the input value of the function, can be scalar or array.
%
%   p : It is the dimension of Von-Mises Fisher distribution which is an
%       integer.
%
% outputs,
%   y : It is the output value of this function, can be scalar or array,
%       depends on the input.
%
% Function is written by Yu-Hui Chen, University of Michigan
    y = besseli(p/2, x)./besseli(p/2-1, x);
end