clc;clear;close all;
M = 8;
K = 8;
H = (randn(K,M)+1i*randn(K,M))*sqrt(0.5);
sigma_q = 1;
Pt = 100;
p = Pt/K*ones(K,1);
sumRate_BC = real(log2(det(eye(K)+H*diag(p)*H')));
W_ZF = ZF(H,Pt);
W_MMSE = MMSE(H,Pt);
SumRate_ZF = myCalSumRate(H,W_ZF);
SumRate_MMSE = myCalSumRate(H,W_MMSE);
Omega = ones(K,1);
iter_max = 200;
eplision = 1e-4;
W0 = randn(M,K)+1i*randn(M,K);
W0 = sqrt(Pt)/norm(W0,'fro')*W0;
[~,SumRate_WSR,~] = WSR(H,W0,Omega,Pt,sigma_q,iter_max,eplision);
function SumRate = myCalSumRate(H,W)
    [K,~] = size(H);
    SumRate = 0;
    for k = 1:K
        hk_H = H(k,:);
        wk = W(:,k);
        numerator = abs(hk_H*wk)^2;
        denominator = sum(norm(hk_H*W)^2)-numerator+1;
        SumRate = SumRate + log2(1+numerator/denominator);
    end
end
