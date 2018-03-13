function [Mu, Kappa, W, logL]=VMFEMz(X1, Pm1,CI1,X2, Pm2,CI2, Num_of_clusters,Num_of_init)
    Pm1 = cat(3, Pm1, -Pm1);
    Pm2 = cat(3, Pm2, -Pm2);
    %%% Duplicate the Euler Angles
    No1 = size(Pm1,3);
    N1 = size(X1,1);
    
    No2 = size(Pm2,3);
    N2 = size(X2,1);
    
    p = size(X1,2);
    total=sum(CI1)+sum(CI2);
    % Precompute the xAp, yAp for invAp_fast
    xAp=0.001:0.1:700;
    yAp = Ap(xAp, p);
    
    % Create container for estimated parameters
    if(nargin<7)
        Num_of_clusters=1;
    end
    if(nargin<8)
        Num_of_init=8;
    end
    
    Mu_All = zeros(Num_of_clusters, p, Num_of_init);
    Kappa_All = zeros(Num_of_clusters, Num_of_init);
    W_All = zeros(Num_of_clusters, Num_of_init);
    R_All = cell(Num_of_init,1);
    L_All = zeros(Num_of_init,1);
    for init=1:Num_of_init   
        % Initialize Mu and Kappa
        R1 = zeros(N1,No1,Num_of_clusters);
        R2 = zeros(N2,No2,Num_of_clusters);
        Mu = zeros(Num_of_clusters, p);
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
                for j=1:No1
                    R1(:,j,clu) = W(clu)*VMFDensity(X1, (Pm1(:,:,j)*Mu(clu,:)')', Kappa(clu));
                end
            end
            
            for clu = 1:Num_of_clusters
                for j=1:No2
                    R2(:,j,clu) = W(clu)*VMFDensity(X2, (Pm2(:,:,j)*Mu(clu,:)')', Kappa(clu));
                end
            end
            % Normalization
            Rdenom = sum(sum(R1,3),2);
            R1 = R1 ./ repmat(Rdenom, [1,No1,Num_of_clusters]);
            
            
            Rdenom = sum(sum(R2,3),2);
            R2 = R2 ./ repmat(Rdenom, [1,No2,Num_of_clusters]);
            %%% M-step
            % estimate W
            W = squeeze(sum(sum(R1,2).*repmat(CI1, [1,1,Num_of_clusters]),1)+sum(sum(R2,2).*repmat(CI2, [1,1,Num_of_clusters]),1));
            W = W / sum(W(:));
            
            for clu = 1:Num_of_clusters
                % estimate Mu
                tmpGamma1 = zeros(No1, p);
                
                for j=1:No1
                    tmpGamma1(j,:) = sum((Pm1(:,:,j)'*X1')'.*repmat(CI1.*R1(:,j,clu), [1 4]));
                end
                Gamma = sum(tmpGamma1,1);
                
                tmpGamma2 = zeros(No2, p);
                for j=1:No2
                    tmpGamma2(j,:) = sum((Pm2(:,:,j)'*X2')'.*repmat(CI2.*R2(:,j,clu), [1 4]));
                end
                Gamma=Gamma+sum(tmpGamma2,1);
                Mu(clu,:) = Gamma / norm(Gamma,2);
                % estimate Kappa
                Kappa(clu) = invAp(norm(Gamma,2)/(W(clu)*total), p, xAp, yAp);
            end
            
            % Calculate the Q function
            Phi1 = zeros(N1,No1,Num_of_clusters);
            for clu = 1:Num_of_clusters
                for j=1:No1
                    Phi1(:,j,clu) = CI1'.*W(clu)*VMFDensity(X1, (Pm1(:,:,j)*Mu(clu,:)')', Kappa(clu));
                end
            end
            Phi2 = zeros(N1,No1,Num_of_clusters);
            for clu = 1:Num_of_clusters
                for j=1:No2
                    Phi2(:,j,clu) = CI2'.*W(clu)*VMFDensity(X2, (Pm2(:,:,j)*Mu(clu,:)')', Kappa(clu));
                end
            end
            L(ite) = sum(log(sum(sum(Phi1+eps,3),2)))+sum(log(sum(sum(Phi2+eps,3),2)));
            
            % update the containers
            Mu_All(:,:,init) = Mu;
            Kappa_All(:,init) = Kappa(:);
            W_All(:,init) = W(:); 
            R_All{init} = R1;
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
        R1 = R_All{dd}; 
        [yy,CIdx] = max(squeeze(sum(R1,2)),[],2);
    else
        CIdx = ones(N1,1);
    end
end