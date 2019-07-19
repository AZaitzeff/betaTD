function [replacemap,betadict,kappa,numgrains]=Germain(neighbors,sizeneigh,dict,beta,IQ,K,tol,minsize)
    max_init=8;
    T=alphatobetatrans();
    Pm=getsymmetries('cubic');
    Pall=zeros(4,4,144);
    for i=1:144
        Pall(:,:,i)=Pm(:,:,ceil(i/6))*T(:,:,mod(i-1,6)+1)';
    end
    betadict=zeros(4,K);
    kappa=zeros(1,K);
    replacemap=zeros(1,K);
    if tol>.75
        tol=tol/180*pi;
    end
    %transfmatrices;
    [~,order]=sort(IQ,'descend');

    %neigh = [[1;0] [-1;0] [0;1] [0;-1]];
    

    z=1;
    for ind=order
        if replacemap(ind)==0
            lenstack=100;
            stack=zeros(1,lenstack);
            stack(1)=ind;
            pop=1;
            next=2;
            while pop<next
               
                cur=stack(pop);
                curori=dict(:,cur);
                curbeta=beta(cur);
                pop=pop+1;
                NB=sizeneigh(cur);
                for nb=1:NB
                    gridind=neighbors(cur,nb);
                    if gridind~=0 && replacemap(gridind)==0 && misorientation(curori,dict(:,gridind),curbeta,beta(gridind))<tol
                        replacemap(gridind)=z;
                        stack(next)=gridind;
                        next=next+1;
                        if next>lenstack
                            templenstack=ceil(lenstack*1.5);
                            tempstack=zeros(1,lenstack);
                            tempstack(1:lenstack)=stack;
                            lenstack=templenstack;
                            stack=tempstack;
                        end
                    end
                end
            end
            if next>minsize
                EBSDtemp=dict(:,stack(1:(next-1)));
                mask=beta(1:(next-1));
                alphaEBSD=EBSDtemp(:,~mask);
                betaEBSD=EBSDtemp(:,mask);
                [mu, kap, ~] = VMFEMzfast(alphaEBSD, Pall,betaEBSD, Pm,max_init,[1;0;0;0],1);
                misoriena=alpbmetric(alphaEBSD,mu);
                misorienb=b2bmetric(betaEBSD,mu);
                med=1.5*median([misoriena misorienb]);
                maska=misoriena<med;
                maskb=misorienb<med;
                [mu, kap, ~] = VMFEMzfast(alphaEBSD(:,maska), Pall,betaEBSD(:,maskb), Pm,1,mu,kap);
                betadict(:,z)=mu;
                kappa(z)=kap;
                for i=1:next-1
                    indtemp=stack(i);
                    if misorientation(mu,dict(:,indtemp),1,beta(indtemp))>tol
                        replacemap(indtemp)=0;
                    end
                end
                z=z+1;
            else
                for i=1:next-1
                    replacemap(stack(i))=0;
                end
            end
            
        end
    end
    
    numgrains=z-1;
    betadict(:,z:K)=[];
    kappa(:,z:K)=[];
end

 function dist=misorientation(x,y,betax,betay)
     if betax==1 && betay==1
         dist=b2bmetric(x,y);
     elseif betax==0 && betay==0
         dist=a2asameparentmetric(x,y);
     elseif betax==0 && betay==1
         dist=alpbmetric(x,y);   
     elseif betax==1 && betay==0
         dist=alpbmetric(y,x);   
     end
 end