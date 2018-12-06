function [Mu, Kappa, logL]=VMFEMzfast(X1, Pm1, X2, Pm2,Num_of_init,Mu,Kappa)
    Pm1 = cat(3, Pm1, -Pm1);
    Pm2 = cat(3, Pm2, -Pm2);
    %%% Duplicate the Euler Angles
    No1 = size(Pm1,3);
    N1 = size(X1,1);
    
    No2 = size(Pm2,3);
    N2 = size(X2,1);
    N=N1+N2;
    p = 4;
    % Precompute the xAp, yAp for invAp_fast
    xAp=0.001:0.1:700;
    yAp = real(Ap(xAp, p));
    
    % Create container for estimated parameters
    if(nargin<5)
        Num_of_init=8;
    end
    if(nargin<6)
        Mu=[];
        Kappa=[];
    end
    
    Mu_All = zeros( Num_of_init,p);
    Kappa_All = zeros(1, Num_of_init);
    L_All = zeros(Num_of_init,1);
    for init=1:Num_of_init   
        % Initialize Mu and Kappa
        R1 = zeros(N1,No1);
        R2 = zeros(N2,No2);
        

        if Num_of_init>1
            Mu = zeros(1, p);
            mu = randn(1,p);
            mu = mu/norm(mu,2);
            Mu(1,:) = mu(:);
            Kappa = 50;
        end
        %a=16;
        %Mu(1,:)=[1,0,0,0];
        %Mu(1,:)=[cos(pi/a),sin(pi/a),0,0];
        %Mu(2,:)=[cos(3*pi/a),sin(3*pi/a),0,0];

        
        
        Num_of_iteration=30;
        L=zeros(Num_of_iteration,1);
        for ite=1:Num_of_iteration

            for j=1:No1
                R1(:,j) = VMFDensityfast(X1, (Pm1(:,:,j)*Mu(1,:)')', Kappa);
            end
            % Normalization
            Rdenom = sum(R1,2);
            R1 = R1 ./ repmat(Rdenom, [1,No1]);

            %%% M-step
            % estimate W

                % estimate Mu
            tmpGamma1 = zeros(No1, p);

            for j=1:No1
                tmpGamma1(j,:) = sum((Pm1(:,:,j)'*X1')'.*repmat(R1(:,j), [1 4]));
            end


            for j=1:No2
                R2(:,j) = VMFDensityfast(X2, (Pm2(:,:,j)*Mu(1,:)')', Kappa);
            end

            % Normalization
            Rdenom = sum(R2,2);
            R2 = R2 ./ repmat(Rdenom, [1,No2]);

            tmpGamma2 = zeros(No2, p);
            for j=1:No2
                tmpGamma2(j,:) = sum((Pm2(:,:,j)'*X2')'.*repmat(R2(:,j), [1 4]));
            end
            Gamma=sum(tmpGamma1,1)+sum(tmpGamma2,1);
            
            Mu(1,:) = Gamma / norm(Gamma,2);
            % estimate Kappa
            Kappa = invAp(norm(Gamma,2)/N, p, xAp, yAp);
            

            Phi1 = zeros(N1,No1);
            for j=1:No1
                Phi1(:,j) = VMFDensityfast(X1, (Pm1(:,:,j)*Mu(1,:)')', Kappa);
            end


            Phi2 = zeros(N2,No2);
            for j=1:No2
                Phi2(:,j) = VMFDensityfast(X2, (Pm2(:,:,j)*Mu(1,:)')', Kappa);
            end
            L(ite)=sum(log(sum(Phi1+eps,2)));+sum(log(sum(Phi2+eps,2)));
            % update the containers
            Mu_All(init,:) = Mu;
            Kappa_All(init) = Kappa;
            L_All(init) = L(ite);
            
            % Check stopping criteria
            if(ite>=2)
                if(abs(L(ite)-L(ite-1))<0.05)
                    break;
                end
            end      
        end
    end
    [~, dd] = max(L_All);
    Mu=Mu_All(dd,:);
    Kappa=Kappa_All(dd);
    logL = L_All(dd);
end