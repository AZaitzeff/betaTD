function [Mu, Kappa, logL]=VMFEMfast(X, Pm,Num_of_initarg,Mus,Kappas)
    Pm = cat(3, Pm, -Pm);
    %%% Duplicate the Euler Angles
    No = size(Pm,3);
    N = size(X,1);
    p = 4;
    
    % Precompute the xAp, yAp for invAp_fast
    xAp=0.001:0.1:700;
    yAp = real(Ap(xAp, p));
    
    % Create container for estimated parameters
    if(nargin<3)
        Num_of_initarg=20;
        Kappas=50;
        Mus=[1,0,0,0];
    end
    Num_of_init=abs(Num_of_initarg);
    Mu_All = zeros( Num_of_init,p);
    Kappa_All = zeros(1, Num_of_init);
    L_All = zeros(Num_of_init,1);
    for init=1:Num_of_init   
        % Initialize Mu and Kappa
        R = zeros(N,No);
        

        if Num_of_initarg>1
            Mu = zeros(1, p);
            mu = randn(1,p);
            mu = mu/norm(mu,2);
            Mu(1,:) = mu(:);
            Kappa = 50;
        else
            Mu = zeros(1, p);
            Mu(1,:) = Mus(init,:);
            Kappa=Kappas;
        end
        %a=16;
        %Mu(1,:)=[1,0,0,0];
        %Mu(1,:)=[cos(pi/a),sin(pi/a),0,0];
        %Mu(2,:)=[cos(3*pi/a),sin(3*pi/a),0,0];

        
        
        Num_of_iteration=30;
        L=zeros(Num_of_iteration,1);
        for ite=1:Num_of_iteration
            %%% E-step
            for j=1:No
                R(:,j) = VMFDensityfast(X, (Pm(:,:,j)*Mu(1,:)')', Kappa);
            end
            % Normalization
            Rdenom = sum(R,2);
            R = R ./ repmat(Rdenom, [1,No]);
            
            %%% M-step
            % estimate W
            
                % estimate Mu
                tmpGamma = zeros(No, p);
                
                for j=1:No
                    tmpGamma(j,:) = sum((Pm(:,:,j)'*X')'.*repmat(R(:,j), [1 4]),1);
                end
                Gamma = sum(tmpGamma,1);
                
                Mu(1,:) = Gamma / norm(Gamma,2);
                % estimate Kappa
                Kappa = invAp(norm(Gamma,2)/N, p, xAp, yAp);
            
            % Calculate the Q function
            Phi = zeros(N,No);
            for j=1:No
                Phi(:,j) = VMFDensityfast(X, (Pm(:,:,j)*Mu(1,:)')', Kappa);
            end
            L(ite) = sum(log(sum(Phi+eps,2)));
            
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
    [yy, dd] = max(L_All);
    Mu=Mu_All(dd,:);
    Kappa=Kappa_All(dd);
    logL = L_All(dd);
end