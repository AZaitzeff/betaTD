function hsym = HexSymmetries()
% HexSymmetries - Quaternions for hexagonal symmetry group.
%
%   USAGE:
%
%   hsym = HexSymmetries
%
%   INPUT:  none
%
%   OUTPUT:
%
%   hsym is 4 x 12,
%        it is the hexagonal symmetry group represented
%        as quaternions
%   
sym= [ 0, 0, 0, 1.;
        1., 0., 0., 0;
        0, 1., 0., 0;
        0, 0, 1., 0;
        sqrt(3)/2,0.,0.,1/2;
        sqrt(3)/2,0.,0.,-1/2;
        1/2,0.,0.,sqrt(3)/2;
        1/2,0.,0.,-sqrt(3)/2;
        0.,1/2,sqrt(3)/2,0.;
        0.,-1/2,sqrt(3)/2,0.;
        0.,sqrt(3)/2,1/2,0.;
        0.,-sqrt(3)/2,1/2,0.;];
%
%  Axis does not need to be normalized in call to QuatOfAngleAxis.
%
hsym = sym';
