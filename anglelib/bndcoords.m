function [xbdcor,ybdcor,sizebdcor,coords,sizecor,minmaxrowcol]=  bndcoords(mapall,K)
[M,N]=size(mapall);


sizecor=ones(K,1);
sizecors=ceil(4*M*N/K);

coords=zeros(K,sizecors);
sizebdcor=ones(K,1);
sizebd=ceil(4*sqrt(M)*sqrt(N)/K);

xbdcor=zeros(K,sizebd);
ybdcor=zeros(K,sizebd);
%coder.varsize('xbdcor','ybdcor','coords');
%bdcor=zeros(1,sizebd);
%neighbors=zeros(K,K);
minmaxrowcol=zeros(K,4);%min row,max row, min col, max col
minmaxrowcol(:,1)=M;
minmaxrowcol(:,2)=0;
minmaxrowcol(:,3)=N;
minmaxrowcol(:,4)=0;

for j=1:N
    for i =1:M
        k=mapall(i,j);
        if k>0
            if sizecor(k)>sizecors
                sizecors=sizecors*2;
                [coords]=growarrayz(coords,sizecor,K,sizecors);
            end
            coords(k,sizecor(k))=i+(j-1)*M;
            sizecor(k)=sizecor(k)+1;
            if i~=1 && i~=M && j~=1 && j~=N
                total=(mapall(i+1,j)==k)+(mapall(i-1,j)==k)+(mapall(i,j+1)==k)+(mapall(i,j-1)==k);
                %neighbors(k,mapall(i+1,j))=1;
                %neighbors(k,mapall(i,j+1))=1;
                if total<4 
                    if sizebdcor(k)>sizebd
                        sizebd=sizebd*2;
                        [xbdcor]=growarrayz(xbdcor,sizebdcor,K,sizebd);
                        [ybdcor]=growarrayz(ybdcor,sizebdcor,K,sizebd);
                        %[bdcor,sizebdcor]=growarray(bdcor,sizebdcor,K,sizebd,-1);
                    end
                    xbdcor(k,sizebdcor(k))=j;
                    ybdcor(k,sizebdcor(k))=i;
                    %bdcor(k,sizebdcor(k))=i+M*(j-1);
                    sizebdcor(k)=sizebdcor(k)+1;
                    if minmaxrowcol(k,1)>i
                        minmaxrowcol(k,1)=i;
                    end
                    if minmaxrowcol(k,2)<i
                        minmaxrowcol(k,2)=i;
                    end
                    if minmaxrowcol(k,3)>j
                        minmaxrowcol(k,3)=j;
                    end
                    if minmaxrowcol(k,4)<j
                        minmaxrowcol(k,4)=j;
                    end
                end
            end
        end
    end
end
sizebdcor(:)=sizebdcor(:)-1;
sizecor(:)=sizecor(:)-1;
%neighbors=(neighbors+neighbors')>.5;
%for k=1:K
%    neighbors(k,k)=0;
%end
end

function [new]=growarrayz(orig,sizeorig,os,ns)
new=zeros(os,ns);
for i=1:os
    v=sizeorig(i)-1;
    new(i,1:v)=orig(i,1:v);
end
end