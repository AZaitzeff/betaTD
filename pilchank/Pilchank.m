function [betadict,done]=Pilchank(neighbors,sizeneigh,dict,betamean,K,vartol,betatol)
    minneigh=2;
    T=alphatobetatrans();
    Pm=getsymmetries('cubic');
    Pall=zeros(4,4,144);
    for i=1:144
        Pall(:,:,i)=Pm(:,:,ceil(i/6))*T(:,:,mod(i-1,6)+1)';
    end
    betadict=zeros(4,K);
    done=zeros(1,K);
    if vartol>.75
        vartol=vartol/180*pi;
    end
    
    if betatol>.75
        betatol=betatol/180*pi;
    end

    for ind=1:K
        ns=sizeneigh(ind);
        if betamean(ind)
            betadict(:,ind)=dict(:,ind);
        else
            if ns>minneigh
                mainalpha=dict(:,ind);
                neighborhood=zeros(1,ns);
                z=1;
                for i=1:ns
                    curind=neighbors(ind,i);
                    if misorientation(mainalpha,dict(:,curind),0,betamean(curind))<vartol
                        neighborhood(z)=curind;
                        z=z+1;
                    end
                end
                neighborhood(z:ns)=[];
                orien=dict(:,neighborhood);
                mask=betamean(neighborhood);
                numalpha=sum(~mask);
                alpha=zeros(4,numalpha+1);
                alpha(:,1)=mainalpha;
                alpha(:,2:numalpha+1)=orien(:,~mask);
                
                beta=orien(:,mask);
                [betasol,unique]=bestbeta(alpha,beta,betatol);
                if unique
                    betadict(:,ind)=betasol(:,1);
                    done(ind)=1;
                    for i=1:(z-1)
                        curind=neighborhood(i);
                        if ~done(curind)
                            betadict(:,curind)=betasol(:,1);
                            done(curind)=1;
                        end
                        
                    end
                end

            end
        end
    end
end

 function dist=misorientation(x,y,betax,betay)
     if betax==1 && betay==1
         dist=b2bmetric(x,y);
     elseif betax==0 && betay==1
         dist=alpbmetric(x,y);   
     elseif betax==1 && betay==0
         dist=alpbmetric(y,x); 
     else
         dist=a2asameparentmetric(x,y);
     end
 end