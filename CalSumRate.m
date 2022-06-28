function SumRate = CalSumRate(H,W,path_loss)
    [K,~] = size(H);
    rho = 10.^(path_loss/10);
    SumRate = 0;
    for k = 1:K
        rhok = rho(k);
        hk_H = H(k,:);
        wk = W(:,k);
        numerator = rhok*abs(hk_H*wk)^2;
        denominator = rhok*sum(norm(hk_H*W)^2)-numerator+1;
        SumRate = SumRate + log2(1+numerator/denominator);
    end
end

