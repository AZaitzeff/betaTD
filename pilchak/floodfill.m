function [map,numgrains]=floodfill(EBSD,beta,IQ,CI,tol,minsize)

    if tol>.75
        tol=tol/180*pi;
    end
    %transfmatrices;
    [m,n,~]=size(EBSD);
    total=m*n;
    map=zeros(m,n);
    valid=CI>1e-4;
    EBSDf=reshape(EBSD, [total 3]);
    quat=E313toq(EBSDf');
    [~,order]=sort(IQ(:),'descend');

    neigh = [[1;0] [-1;0] [0;1] [0;-1] [1;1] [-1;1] [-1;-1] [1;-1] [2;0] [-2;0] [0;2] [0;-2]];
    ind2sub1 = @(k)[rem(k-1, m)+1; (k - rem(k-1, m) - 1)/m + 1];
    sub2ind1 = @(u)(u(2)-1)*m+u(1);
    %neigh = [[1;0] [-1;0] [0;1] [0;-1]];
    [~,NB]=size(neigh);
    
    N=numel(order);

    z=1;
    for i=1:N
        ind=order(i);
        if map(ind)==0 && valid(ind)
            map(ind)=z;
            lenstack=100;
            stack=zeros(1,lenstack);
            %coder.varsize('stack',[1 500^2],[0,1]); %Uncomment if have coder
            stack(1)=ind;
            pop=1;
            next=2;
            while pop<next
               
                cur=stack(pop);
                curori=quat(:,cur);
                curbeta=beta(cur);
                pop=pop+1;
                for nb=1:NB
                    gridind=findneighbor(cur,nb);
                    if gridind~=0 && map(gridind)==0 && valid(ind) && misorientation(curori,quat(:,gridind),curbeta,beta(gridind))<tol
                        map(gridind)=z;
                        stack(next)=gridind;
                        next=next+1;
                        if next>lenstack
                            templenstack=ceil(lenstack*1.5);
                            tempstack=zeros(1,templenstack);
                            tempstack(1:lenstack)=stack;
                            lenstack=templenstack;
                            stack=tempstack;
                        end
                    end
                end
            end
            if next>minsize
                z=z+1;
            else
                for i=1:next-1
                    map(stack(i))=0;
                end
            end
            
        end
    end
    numgrains=z-1;
    function ind=findneighbor(k,i)
    vec=ind2sub1(k)+neigh(:,i);
    if vec(1)<=m && vec(1)>=1 && vec(2)<=n && vec(2)>=1
        ind=sub2ind1( ind2sub1(k)+neigh(:,i) );
    else
        ind=0;
    end
end
end


 function dist=misorientation(x,y,betax,betay)
     if betax==1 && betay==1
         dist=b2bmetric(x,y);
     elseif betax==0 && betay==0
         dist=a2ametric(x,y);
     else
         dist=inf;
     end
 end
