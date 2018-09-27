function [bMu,bxc,br]=kmeansEBSD(EBSD,CI,x,lambda,K,num_iter,MAXITER)
        Pm=getsymmetries('cubic');
        numsub=200;
        [N,z]=size(EBSD);
        if nargin<6
            num_iter=8;
        end
        if nargin<7
            MAXITER=20;
        end

        bval=inf;
        for iter =1:num_iter
        Mu=zeros(K,z);
        Kappas=zeros(K,1);
        xc=rand(2,K);
        xc(1,:)=xc(1,:).*(max(x(:,1))-min(x(:,1)))+min(x(:,1));
        xc(2,:)=xc(2,:).*(max(x(:,2))-min(x(:,2)))+min(x(:,2));
        
        rl=zeros(N,1);
        val=zeros(N,K);
        for k=1:K
            mean=sum((x-xc(:,k)').^2,2);
            val(:,k)=mean;
        end
        [~,r]=min(val,[],2);
        for k=1:K
            mask=(r==k);
            xc(1,k)=sum(x(mask,1))/sum(mask);
            xc(2,k)=sum(x(mask,2))/sum(mask);
            newind=datasample(find(mask),numsub,'Weights',CI(mask));
            CItemp=ones(size(CI(newind)));
            [mu, kappa, ~, ~,~] = VMFEM(EBSD(newind,:), Pm,CItemp,1);
            Mu(k,:)=mu;
            Kappas(k) = kappa;
        end
        
        for t=1:MAXITER
            for k=1:K
                mean=sum((x-xc(:,k)').^2,2);
                rot=alpbmetric(EBSD,Mu(k,:));
                val(:,k)=rot'+lambda*mean;
                
            end
            [~,r]=min(val,[],2);

            for k=1:K
                mask=(r==k);
                xc(1,k)=sum(x(mask,1))/sum(mask);
                xc(2,k)=sum(x(mask,2))/sum(mask);
                newind=datasample(find(mask),numsub,'Weights',CI(mask));
                CItemp=ones(size(CI(newind)));
                [mu, kappa, ~, ~,~] = VMFEM(EBSD(newind,:), Pm,CItemp,1,1,Mu(k,:),Kappas(k));
                Mu(k,:)=mu;
                Kappas(k) = kappa;
            end
            if(sum(abs(rl-r))<.05)
                break
            end
            rl=r;

        end
            val=0;
            for k=1:K
                mask=(r==k);
                mean=sum((x(mask,:)-xc(:,k)').^2,2);
                rot=alpbmetric(EBSD(mask,:),Mu(k,:));
                val=val+sum(rot'+lambda*mean);
                
            end
            if val<bval
                bval=val;
                bMu=Mu;
                bxc=xc;
                br=r;
            end
        end