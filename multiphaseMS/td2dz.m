function [levelsets] = td2dz(levelsets,dt,flag)
if nargin<3
    flag=0;
end
h=1/100;
N = size(levelsets,2); % Number of phases.
[n1,n2] = size(levelsets{1});
if flag
    nx=2*n1;
    ny=2*n2;
else
    border1=10;
    border2=30;
    border=border1+border2;
    nx=n1+2*border;
    ny=n2+2*border;
    
end

%Z = -ones(n1,n2); % Aux. variable, used in calling ls2vf2D.
%ind = find(Z<0);  % Aux. variable, used in calling ls2vf2D.
%[sub1,sub2] = ind2sub([n1 n2],ind); % Aux. variable, used in calling ls2vf2D.
%indicator = 0*Z;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nu = nx;
nv = ny;
du = 1/(nu*h);
dv = 1/(nv*h);
u = [-nu/2:nu/2-1]*du;
v = [-nv/2:nv/2-1]*dv;
[uu,vv] = ndgrid(u,v);
K=exp(-(uu.^2 + vv.^2)*(4*pi^2*dt));
K=ifftshift(K);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Matrix of surface tensions:
%S = ones(N,N) - eye(N,N); % Equal surface tensions.

% Time stepping.

  % Comparison functions initialized:


  % Convolutions:
  for k=1:N
    %indicator = 0*Z;
    %vals = ls2vf2D(int32(sub1),int32(sub2),levelsets{k},Z);
    %indicator(ind) = vals;
    indicator=levelsets{k}>0;
    if flag
        indicator=padarray(indicator,[n1 n2],'symmetric','post');
    else
        indicator=padarray(indicator,[border1 border1],'replicate');
        indicator=padarray(indicator,[border2 border2]);
    end
    
    convolution = real(ifft2( fft2(indicator) .* K ));
    % Convolution dropped into comparison functions it appears in:
    if flag
        levelsets{k} =convolution(1:n1,1:n2);
    else
        levelsets{k} =convolution(border+1:end-border,border+1:end-border);
    end
  end
 % for t. Time stepping ends.