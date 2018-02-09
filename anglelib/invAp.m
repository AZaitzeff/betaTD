function [x] = invAp_fast(y, p, xAp, yAp)
%   This function is used for calculate the density of Von-Mises Fisher
%   distribution. It is basically the Inverse function of Ap and the value
%   is calculated by interpolating the pre-computed Ap function value.
% 
% inputs,
%   y : It is the input value of the function, can be scalar or array.
%
%   p : It is the dimension of Von-Mises Fisher distribution which is an
%       integer.
%
% outputs,
%   x : It is the output value of this function, can be scalar or array,
%       depends on the input.
%
% Function is written by Yu-Hui Chen, University of Michigan

% Find the closest one 
x = interp1(yAp, xAp, y, 'linear', 'extrap');
    
end