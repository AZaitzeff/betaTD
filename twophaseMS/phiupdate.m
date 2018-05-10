function [uo,g1,g2] = phiupdate(nt,dt,uin,zfun,cnstrt,EBSD,CI,betas,fid,g1,g2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [uo] = phiupdatecombined(nt,dt,uin,fEBSD,fBSE,fidEBSD,fidBSE,radius)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Piecewise constant Mumford-Shah model with filter (circle) based on Chan and Vese.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input:

%   nt = Number of time steps to take.

%   dt = Time step size.

%   uin = Initial level set function.

%   fEBSD = Blurred image to be segmented (EBSD data).

%   fBSE = sharp image to be used (BSE data)

%   fidEBSD = Fidelity constant for EBSD data.

%   fidBSE =  Fidelity constant for BSE data.

%   radius = radius of circle filter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Output:

%   uo = Final level set function.

%   g1 = Best unknown constant Kikuchi pattern 1 for uo. (EBSD)

%   g2 = Best unknown constant Kikuchi pattern 2 for uo. (EBSD)

%   c1 = Best unknown constant 1 for uo. (BSE)

%   c2 = Best unknown constant for uo. (BSE)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

me = 0.1;  % Regularization of denominators.
T=alphatobetatrans();
Pm=getsymmetries('cubic');
Pall=zeros(4,4,144);
for i=1:144
    Pall(:,:,i)=T(:,:,mod(i-1,6)+1)'*Pm(:,:,ceil(i/6));
end
[m,n] = size(uin); %size of level set function
[~,~,z]=size(EBSD);
EBSDflat=reshape(EBSD,[m*n,z]);
CIflat=reshape(CI,[m*n,1]);
CIb = buffer2(CI);
EBSDflat=E313toq(EBSDflat);
betasflat=reshape(betas,[m*n,1]);
betamask=betasflat==1;
alphamask=betasflat==0;
invh=100;
tol=.01; %measuring distance
area=sum(zfun(:));
u = buffer2(uin);   % Size is now (m+4)x(n+4). adds zero derivative boudary conditions
zfunb = padarray(zfun,[2,2]);
cnstrtb = buffer2(cnstrt);
interval=2000;
region0=zfun.*(u(3:m+2,3:n+2)>.5);
region1=zfun.*(u(3:m+2,3:n+2)>.5);
newregion1=region1;
flag=0;
numsub=500;
newg1=g1;
newg2=g2;

for t = 1:nt % Main time loop

    % Difference quotients:  

    u_x = forwardx(u)*invh;

    u_y = forwardy(u)*invh;  

    Du = sqrt( u_x.^2 + u_y.^2 + me );
    if mod(t,interval)==1
        
        mask1=zfun.*(u(3:m+2,3:n+2)>.5);
        mask2=zfun.*(u(3:m+2,3:n+2)<=.5);
        
        check1=sum(mask1(:));
        check2=sum(mask2(:));
        
        if check1>7 && check2>7
            %Replace
            val=sum(sum(abs(newregion1-region1)))/area;
            if val>.05
                indices=find(mask1(:));

                newind=datasample(indices,numsub,'Weights',CIflat(indices));
                CItemp=ones(size(CIflat(newind)));
                EBSDtemp=EBSDflat(newind,:);
                maskalpha=alphamask(newind);
                maskbeta=betamask(newind);
                if sum(maskbeta)==0
                    [newg1, ~, ~, ~] = VMFEM(EBSDtemp, Pall,CItemp);
                elseif sum(maskalpha)==0
                    [newg1, ~, ~, ~] = VMFEM(EBSDtemp, Pm,CItemp);
                else
                    [newg1, ~, ~, ~] = VMFEMz(EBSDtemp(maskalpha,:), Pall,CItemp(maskalpha),...
                        EBSDtemp(maskbeta,:), Pm,CItemp(maskbeta));
                end
                indices=find(mask2(:));

                newind=datasample(indices,numsub,'Weights',CIflat(indices));
                CItemp=ones(size(CIflat(newind)));
                EBSDtemp=EBSDflat(newind,:);
                maskalpha=alphamask(newind);
                maskbeta=betamask(newind);
                if sum(maskbeta)==0
                    [newg2, ~, ~, ~,~] = VMFEM(EBSDtemp, Pall,CItemp);
                elseif sum(maskalpha)==0
                    [newg2, ~, ~, ~,~] = VMFEM(EBSDtemp, Pm,CItemp);
                else
                    [newg2, ~, ~, ~] = VMFEMz(EBSDtemp(maskalpha,:), Pall,CItemp(maskalpha),...
                        EBSDtemp(maskbeta,:), Pm,CItemp(maskbeta));
                end
                
                region1=zfun.*(u(3:m+2,3:n+2)>.5);
            end
        elseif t<(interval*1.5)
            indices=find(zfun(:));

            newind=datasample(indices,numsub,'Weights',CIflat(indices));
            CItemp=ones(size(CIflat(newind)));
            EBSDtemp=EBSDflat(newind,:);
            maskalpha=alphamask(newind);
            maskbeta=betamask(newind);
            if sum(maskbeta)==0
                [mu, ~, ~, ~,~] = VMFEM(EBSDtemp, Pall,CItemp,2,20);
            elseif sum(maskalpha)==0
                [mu, ~, ~, ~,~] = VMFEM(EBSDtemp, Pm,CItemp,2,20);
            else
                [mu, ~, ~, ~] = VMFEMz(EBSDtemp(maskalpha,:), Pall,CItemp(maskalpha),...
                    EBSDtemp(maskbeta,:), Pm,CItemp(maskbeta),2,20);
            end
             
            newg1=mu(1,:);
            newg2=mu(2,:);
            newg1=newg1/norm(newg1);
            newg2=newg2/norm(newg2);
            if sum(maskbeta)==0
                val1=sum(CItemp'.*alpbmetric(EBSDtemp,newg1));
                val2=sum(CItemp'.*alpbmetric(EBSDtemp,newg2));

            elseif sum(maskalpha)==0
                val1=sum(CItemp'.*b2bmetric(EBSDtemp,newg1));
                val2=sum(CItemp'.*b2bmetric(EBSDtemp,newg2));
            else
                val1=sum(CItemp(maskalpha)'.*alpbmetric(EBSDtemp(maskalpha,:),newg1));
                val2=sum(CItemp(maskalpha)'.*alpbmetric(EBSDtemp(maskalpha,:),newg2));
                val1=val1+sum(CItemp(maskbeta)'.*b2bmetric(EBSDtemp(maskbeta,:),newg1));
                val2=val2+sum(CItemp(maskbeta)'.*b2bmetric(EBSDtemp(maskbeta,:),newg2));
            end
                
            if val1>val2
                temp=newg1;
                newg1=newg2;
                newg2=temp;
            end
            %ratio=(sum(CIflat(newind(Cind==1))'.*alpbmetric(EBSDflat(newind(Cind==1),:),newg1))+...
            %    sum(CIflat(newind(Cind==2))'.*alpbmetric(EBSDflat(newind(Cind==2),:),newg2)))/...
            %    sum(CIflat(newind)'.*alpbmetric(EBSDflat(newind,:),g1));
            %if ratio>.85
            %    break
            %end
            
            %newg1=newg1';
            %newg2=newg2';

        end
        if t==1||b2bmetric(g1,newg1)>tol || b2bmetric(g2,newg2)>tol
            g1=newg1;
            g2=newg2;
            X=zeros(m,n);
            if sum(betamask)>0
                X(betamask)=b2bmetric(EBSDflat(betamask,:),g1)-b2bmetric(EBSDflat(betamask,:),g2);
            end
            if sum(alphamask)>0
                X(alphamask)=alpbmetric(EBSDflat(alphamask,:),g1)-alpbmetric(EBSDflat(alphamask,:),g2);
            end
            Xb=buffer2(X);
            region0=zfun.*(u(3:m+2,3:n+2)>.5);
            flag=0;
        else
            newregion=zfun.*(u(3:m+2,3:n+2)>.5);
            if sum(sum(abs(newregion-region0)))==0
                if flag
                    break
                else
                    flag=1;
                end
            else
                region0=newregion;
                flag=0;
            end
        end

    end
    


    rhs = (backwardx(zfunb.*u_x./Du)*invh + backwardy(zfunb.*u_y./Du)*invh)-fid*CIb.*Xb;
    
    u = u + dt*rhs;
    u(u<0)=0;
    u(u>1)=1;
    u=u.*cnstrtb;
    u=upbuffer2(u);
    newregion1=zfun.*(u(3:m+2,3:n+2)>.5);
end



uo = u(3:m+2,3:n+2);


function [ubo]=upbuffer2(ub)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function updates the buffer layer of a buffered matrix, 

% using extension, where the buffer layer thickness is 2. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Notes: No error checking is built in. Make sure size of input

%    matrix is large enough: both dims. > 4 needed.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



m=size(ub,1)-4;

n=size(ub,2)-4;



% Update the buffer:



%ub(3:m+2,1:2)=ub(3:m+2,n+1:n+2);

%ub(3:m+2,n+3:n+4)=ub(3:m+2,3:4);

%ub(1:2,1:n+4)=ub(m+1:m+2,1:n+4);

%ub(m+3:m+4,1:n+4)=ub(3:4,1:n+4);

ub(3:m+2,1)=ub(3:m+2,3);
ub(3:m+2,2)=ub(3:m+2,3);
ub(3:m+2,n+3)=ub(3:m+2,n+2);
ub(3:m+2,n+4)=ub(3:m+2,n+2);
ub(1,1:n+4)=ub(3,1:n+4);
ub(2,1:n+4)=ub(3,1:n+4);
ub(m+3,1:n+4)=ub(m+2,1:n+4);
ub(m+4,1:n+4)=ub(m+2,1:n+4);

ubo=ub;



function [ub]=buffer2(u)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function buffers up a given matrix with a buffer layer

% of thickness 2. It uses periodic boundary conditions to fill

% the buffer layer.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Uses: upbuffer2.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Notes: If the input matrix u has size mxn, then the output 

%    matrix ub has size (m+4)x(n+4).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



[m,n]=size(u);



ub=zeros(m+4,n+4);  % Initialize ub.

ub(3:m+2,3:n+2)=u;  % Bulk part of ub is just u.



ub=upbuffer2(ub);   % Fills in buffer layer via extension.



function [dy]=backwardy(u)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Backward difference operator in the y-direction for the input

% matrix u, which is assumed to be buffered (layer thickness 2).

% Periodic boundary conditions used.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Uses: upbuffer2.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



[m,n]=size(u);  % Buffered size.

m=m-4;      % Find the unbuffered size.

n=n-4;      % Find the unbuffered size.



dy=u;       % Initialization for dy.



% Carry out the operation on unbuffered part of u: 

dy(3:m+2,3:n+2)=( u(3:m+2,3:n+2) - u(2:m+1,3:n+2) );



% Update the buffer region:

dy=upbuffer2(dy);



function [dx]=backwardx(u)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Backward difference operator in the x-direction for the input

% matrix u, which is assumed to be buffered (layer thickness 2).

% Periodic boundary conditions used.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Uses: upbuffer2.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



[m,n]=size(u);  % Buffered size.

m=m-4;      % Find the unbuffered size.

n=n-4;      % Find the unbuffered size.



dx=u;       % Initialization for dx.



% Carry out operation on unbuffered part of u:

dx(3:m+2,3:n+2)=( u(3:m+2,3:n+2) - u(3:m+2,2:n+1) );



% Update the buffer region:

dx=upbuffer2(dx);





function [dy]=forwardy(u)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Forward difference operator in the y-direction for the input

% matrix u, which is assumed to be buffered (layer thickness 2).

% Periodic boundary conditions used.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Uses: upbuffer2.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



[m,n]=size(u);  % Buffered size.

m=m-4;      % Find the unbuffered size.

n=n-4;      % Find the unbuffered size.



dy=u;       % Initialization for dy.



% Carry out operation on unbuffered part of u:

dy(3:m+2,3:n+2)=( u(4:m+3,3:n+2) - u(3:m+2,3:n+2) );



% Update the buffer region:

dy=upbuffer2(dy);





function [dx]=forwardx(u)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Forward difference operator in the x-direction for the input 

% matrix u, which is assumed to be buffered (layer thickness 2). 

% Periodic boundary conditions used.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Uses: upbuffer2.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



[m,n]=size(u);  % Buffered size.

m=m-4;      % Find the unbuffered size.

n=n-4;      % Find the unbuffered size.



dx=u;       % Initialization for dx.



% Carry out operation on unbuffered part of u:

dx(3:m+2,3:n+2)=( u(3:m+2,4:n+3) - u(3:m+2,3:n+2) );



% Update the buffer region:

dx=upbuffer2(dx);



function y = DH(x,me)



d = 0.0001;

y = (H(x+d,me)-H(x-d,me))/(2*d);



function y = H(x,me)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function y = H(x,me)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Approximate Heaviside function.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



y = (tanh(x/me)+1)/2;