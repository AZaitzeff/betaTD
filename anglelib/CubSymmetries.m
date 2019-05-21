function csym = CubSymmetries()
% CubSymmetries - Return quaternions for cubic symmetry group.
%
%   USAGE:
%
%   csym = CubSymmetries
%
%   INPUT:  none
%
%   OUTPUT:
%
%   csym is 4 x 24, 
%        quaternions for the cubic symmetry group
%
sym= [ 0, 0, 0, 1.;
        1., 0., 0., 0;
        0, 1., 0., 0;
        0, 0, 1., 0;
        sqrt(2)/2,sqrt(2)/2, 0, 0;
        sqrt(2)/2,0, sqrt(2)/2, 0;
        sqrt(2)/2,0, 0, sqrt(2)/2;
        sqrt(2)/2,-sqrt(2)/2, 0,  0;
        sqrt(2)/2,0, -sqrt(2)/2, 0;
        sqrt(2)/2,0, 0, -sqrt(2)/2;
        0,sqrt(2)/2, sqrt(2)/2, 0;
        0,-sqrt(2)/2, sqrt(2)/2, 0;
        0,0, sqrt(2)/2, sqrt(2)/2;
        0,0, sqrt(2)/2, -sqrt(2)/2;
        0,sqrt(2)/2, 0, sqrt(2)/2;
        0,sqrt(2)/2, 0, -sqrt(2)/2;
        1/2, 1/2, 1/2, 1/2;
        1/2, -1/2, -1/2, -1/2;
        1/2, -1/2, 1/2, 1/2;
        1/2, 1/2, -1/2, 1/2;
        1/2, 1/2, 1/2, -1/2;
        1/2, 1/2, -1/2, -1/2;
        1/2, -1/2, 1/2, -1/2;
        1/2, -1/2, -1/2,1/2;];
%
%  Axis does not need to be normalized; it is done
%  in call to QuatOfAngleAxis.
%
csym = sym';
