function [Mu, Kappa, W,xc,xs, logL, CIdx]=VMFEMwX(X, Pm,CI,lambda,pos, Num_of_clusters,Num_of_init)
    Pm = cat(3, Pm, -Pm);
    %%% Duplicate the Euler Angles
    No = size(Pm,3);
    N = size(X,1);
    p = size(X,2);
    total=sum(CI);
    % Precompute the xAp, yAp for invAp_fast
    xAp=0.001:0.1:700;
    yAp = Ap(xAp, p);
    
    % Create container for estimated parameters
    if(nargin<6)
        Num_of_clusters=1;
    end
    if(nargin<7)
        Num_of_init=8;
    end
    
    Mu_All = zeros(Num_of_clusters, p, Num_of_init);
    Kappa_All = zeros(Num_of_clusters, Num_of_init);
    W_All = zeros(Num_of_clusters, Num_of_init);
    R_All = cell(Num_of_init,1);
    L_All = zeros(Num_of_init,1);
    for init=1:Num_of_init   
        % Initialize Mu and Kappa
        R = zeros(N,No,Num_of_clusters);
        Mu = zeros(Num_of_clusters, p);
        xc=rand(2,Num_of_clusters);
        xc(1,:)=xc(1,:).*(max(pos(:,1))-min(pos(:,1)))+min(pos(:,1));
        xc(2,:)=xc(2,:).*(max(pos(:,2))-min(pos(:,2)))+min(pos(:,2));
        xs=zeros(2,2,Num_of_clusters);
        xs(1,1,:)=1;
        xs(2,2,:)=1;
        for clu = 1:Num_of_clusters
            mu = randn(1,p);
            mu = mu/norm(mu,2);
            Mu(clu,:) = mu(:)';
        end
        %a=16;
        %Mu(1,:)=[1,0,0,0];
        %Mu(1,:)=[cos(pi/a),sin(pi/a),0,0];
        %Mu(2,:)=[cos(3*pi/a),sin(3*pi/a),0,0];

        Kappa = 50*ones(Num_of_clusters,1);
        W = (1/Num_of_clusters)*ones(Num_of_clusters,1);
        
        Num_of_iteration=30;
        L=zeros(Num_of_iteration,1);
        for ite=1:Num_of_iteration
            %%% E-step
            for clu = 1:Num_of_clusters
                mean=pos-xc(:,clu)';
                temp=1/(2*pi*sqrt(det(squeeze(xs(:,:,clu)))))*exp(-1/2*sum(mean.*(squeeze(xs(:,:,clu))\mean')',2));
                for j=1:No
                    R(:,j,clu) = W(clu)*VMFDensity(X, (Pm(:,:,j)*Mu(clu,:)')', Kappa(clu))+lambda/No*temp;
                end
            end
            % Normalization
            Rdenom = sum(sum(R,3),2);
            R = R ./ repmat(Rdenom, [1,No,Num_of_clusters]);
            
            %%% M-step
            % estimate W
            r=squeeze(sum(R,2));
            rdenom=sum(r,2);
            
            r = r./repmat(rdenom, [1 Num_of_clusters]);
            
            W = squeeze(sum(sum(R,2).*repmat(CI, [1,1,Num_of_clusters]),1));
            
            W = W / sum(W(:));
            
            
            
            for clu = 1:Num_of_clusters
                % estimate Mu
                tmpGamma = zeros(No, p);
                
                for j=1:No
                    tmpGamma(j,:) = sum((Pm(:,:,j)'*X')'.*repmat(CI.*R(:,j,clu), [1 4]));
                end
                Gamma = sum(tmpGamma,1);
                
                Mu(clu,:) = Gamma / norm(Gamma,2);
                % estimate Kappa
                Kappa(clu) = invAp(norm(Gamma,2)/(W(clu)*total), p, xAp, yAp);
                
                temp=zeros(size(pos));
                temp(:,1)=CI.*r(:,clu).*pos(:,1);
                temp(:,2)=CI.*r(:,clu).*pos(:,2);
                xc(1,clu)=sum(temp(:,1))/sum(CI.*r(:,clu));
                xc(2,clu)=sum(temp(:,2))/sum(CI.*r(:,clu));
                sigma=zeros(2,2);
                for i=1:N
                    sigma=sigma+temp(i,:)'*pos(i,:);
                end
                xs(:,:,clu)=sigma/sum(CI.*r(:,clu))-xc(:,clu)*xc(:,clu)';
            end
            
            % Calculate the Q function
            Phi = zeros(N,No,Num_of_clusters);
            for clu = 1:Num_of_clusters
                %size(r(:,clu))
                mean=pos-xc(:,clu)';
                %size(sum(mean.*(squeeze(xs(:,:,clu))\mean')',2))
                temp=W(clu)/(2*pi*sqrt(det(squeeze(xs(:,:,clu)))))*r(:,clu).*exp(-1/2*sum(mean.*(squeeze(xs(:,:,clu))\mean')',2));
                for j=1:No
                    Phi(:,j,clu) = W(clu)*CI.*VMFDensity(X, (Pm(:,:,j)*Mu(clu,:)')', Kappa(clu))+lambda/No*CI.*temp;
                end
            end
            L(ite) = sum(log(sum(sum(Phi+eps,3),2)));
            
            % update the containers
            Mu_All(:,:,init) = Mu;
            Kappa_All(:,init) = Kappa(:);
            W_All(:,init) = W(:); 
            R_All{init} = R;
            L_All(init) = L(ite);
            
            % Check stopping criteria
            if(ite>=2)
                if(abs(L(ite)-L(ite-1))<0.05)
                    break;
                end
            end
        end       
    end
    
    [yy, dd] = max(L_All);
    Mu=Mu_All(:,:,dd);
    Kappa=Kappa_All(:,dd);
    W = W_All(:,dd);
    logL = L_All(dd);
    if(Num_of_clusters>1)
        R = R_All{dd}; 
        [yy,CIdx] = max(squeeze(sum(R,2)),[],2);
    else
        CIdx = ones(N,1);
    end
end