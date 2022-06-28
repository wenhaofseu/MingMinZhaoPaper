function [W,rate,rate_all] = WSR(H,W,Omega,Pt,sigma_q,iter_max,eplision)
    %W = [w1,w2,...,wK];
    %H = [h_1^H;h_2^H;...;h_K^H]
    [K,M] = size(H);
    chi = zeros(K,1);
    kappa = zeros(K,1);
    rate_all = zeros(iter_max,1);
    for iter = 1:iter_max
        for k = 1:K
            hk_H = H(k,:);
            wk = W(:,k);
            %% update chi
            chi(k) = 1/(norm(hk_H*W)^2+sigma_q)*hk_H*wk;
            %% update kappa
            kappa(k) = 1/(1-conj(chi(k))*hk_H*wk);
        end
        temp1 = zeros(M);
        temp2 = zeros(M);
        for i = 1:K
            hi_H = H(i,:);
            temp1 = temp1 + Omega(i)*abs(chi(i))^2*kappa(i)*(hi_H'*hi_H);
            Omega(i)^2*abs(chi(i))^2*abs(kappa(i))^2*(hi_H'*hi_H);
            temp2 = temp2 + Omega(i)^2*abs(chi(i))^2*abs(kappa(i))^2*(hi_H'*hi_H);
        end
        [D,V]=eig(temp1);%D*V*D' = temp1
        PHI = D'*temp2*D;
        %% check lambda = 0 is optimal?
        flag = 0;
        if rank(temp1)==M
            temp3 = 0;
            for m = 1:M
                temp3 = temp3+PHI(m,m)/V(m,m)^2;
            end
            if temp3<=Pt
                flag = 1;
            end
        end
        %% calculate lambda
        lambda_min = 0;
        lambda_max = sqrt(trace(PHI)/Pt);
        while(flag==0)
            lambda = (lambda_min+lambda_max)/2;
            temp4 = 0;
            for m = 1:M
                temp4 = temp4+PHI(m,m)/(V(m,m)+lambda)^2;
            end
            temp4 = real(temp4);
            if(temp4<Pt)
                lambda_max = lambda;
            else
                lambda_min = lambda;
            end
            %% bisection converage
            if abs(lambda_max-lambda_min)<1e-3 && abs(temp4-Pt)<1e-5
                break
            elseif abs(lambda_max-lambda_min)<1e-7 && lambda_min==0
                break
            end
        end
        %% update wk
        for k = 1:K
            hk_H = H(k,:);
            wk = Omega(k)*chi(k)*kappa(k)*((temp1+lambda*eye(M))\hk_H');
            W(:,k) = wk;
        end
        norm(W,'fro')^2
        %% calculate rate
        rate = 0;
        for k = 1:K
            hk_H = H(k,:);
            wk = W(:,k);
            gamma_k = abs(hk_H*wk)^2/(norm(hk_H*W)^2-abs(hk_H*wk)^2+sigma_q);
            rate = rate + Omega(k)*log2(1+gamma_k);
        end
        rate_all(iter) = rate;
        %% algorithm converage
        if iter>=5 && abs(rate-rate_all(iter-1))/abs(rate)<eplision
            rate_all = rate_all(1:iter);
            break
        end
    end
end

