clc;clear;close all;
M = 32;
K = 8;
rng(1);
H = sqrt(1/2)*(randn(K,M)+1i*randn(K,M));
Omega = ones(K,1);
Pt = 10;
sigma_q = 1;
iter_max = 200;
eplision = 1e-4;
for num = 1:100
    num
    W = randn(M,K)+1i*randn(M,K);
    W = sqrt(Pt)/norm(W,'fro')*W;
    [W,rate,rate_all] = WSR(H,W,Omega,Pt,sigma_q,iter_max,eplision);
    plot(1:length(rate_all),rate_all);hold on;
    xlabel('iteration number');
    ylabel('Sum Rate');
end
%% calculate rate
W_ZF=H'/(H*H');
W_ZF = sqrt(Pt)/norm(W_ZF,'fro')*W_ZF;
W_MMSE = H'/(H*H'+K*sigma_q /Pt*eye(K));
W_MMSE = sqrt(Pt)/norm(W_MMSE,'fro')*W_MMSE;
rate_ZF = 0;
rate_MMSE= 0;
for k = 1:K
    hk_H = H(k,:);
    wk_ZF = W_ZF(:,k);
    gamma_k_ZF = abs(hk_H*wk_ZF)^2/(norm(hk_H*W_ZF)^2-abs(hk_H*wk_ZF)^2+sigma_q);
    wk_MMSE = W_MMSE(:,k);
    gamma_k_MMSE = abs(hk_H*wk_MMSE)^2/(norm(hk_H*W_MMSE)^2-abs(hk_H*wk_MMSE)^2+sigma_q);
    rate_ZF = rate_ZF + Omega(k)*log2(1+gamma_k_ZF);
    rate_MMSE = rate_MMSE + Omega(k)*log2(1+gamma_k_MMSE);
end