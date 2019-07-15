function [colors,bnds]=  bndmap(mapall,dict,thres)
%{
Used for plotting
%}
K=size(dict,1);
[M,N]=size(mapall);
colors=ones(M,N,3);
bnds=zeros(M,N);
map=zeros(K,K);
thres=thres/180*pi;
if numel(thres)==0
for j=1:N
    for i =1:M
        k=mapall(i,j);
        if i~=M && mapall(i+1,j)~=k
            bnds(i,j)=1;
            colors(i,j,:)=0;
        elseif i~=1 && mapall(i-1,j)~=k
            bnds(i,j)=1;
            colors(i,j,:)=0;
        elseif j~=N && mapall(i,j+1)~=k
            bnds(i,j)=1;
            colors(i,j,:)=0;
        elseif j~=1 && mapall(i,j-1)~=k
            bnds(i,j)=1;
            colors(i,j,:)=0;
        end

    end
end
    
else
for k1=1:K
    for k2=(k1+1):K
        temp=b2bmetric(dict(k1,:),dict(k2,:));
        map(k1,k2)=temp;
        map(k2,k1)=temp;
    end
end
colormapz=[[0,0,0];[0.9290, 0.6940, 0.1250];[0.8500, 0.3250, 0.0980];...
    [0.4660, 0.6740, 0.1880];[0.3010, 0.7450, 0.9330];[0.6350, 0.0780, 0.1840];[0.4940, 0.1840, 0.5560]];
for j=1:N
    for i =1:M
        k=mapall(i,j);
        if i~=M && mapall(i+1,j)~=k
            bnds(i,j)=1;
            colors(i,j,:)=colorfind(map(k,mapall(i+1,j)),colormapz,thres);
        elseif i~=1 && mapall(i-1,j)~=k
            bnds(i,j)=1;
            colors(i,j,:)=colorfind(map(k,mapall(i-1,j)),colormapz,thres);
        elseif j~=N && mapall(i,j+1)~=k
            bnds(i,j)=1;
            colors(i,j,:)=colorfind(map(k,mapall(i,j+1)),colormapz,thres);
        elseif j~=1 && mapall(i,j-1)~=k
            bnds(i,j)=1;
            colors(i,j,:)=colorfind(map(k,mapall(i,j-1)),colormapz,thres);
        end

    end
end
end
end

function color=colorfind(map,colors,thres)
    Z=numel(thres);
    cind=1;
    for z=1:Z
        if map<thres(z)
            cind=(z+1);
            break
        end
    end
    color=colors(cind,:);

end
