clc;clear;close all;
number = 1000;%number of Monter Carlo
M = 8;
K = 8;
P_dB_range = -10:5:20;
Prange = 10.^(P_dB_range/10);
path_loss = 20+unifrnd(-5,5,1,K);%noise power:-80dBm;average path loss:130dB
rho = 10.^(path_loss/10);
SumRate_ZF = zeros(1,length(Prange));
SumRate_MMSE = zeros(1,length(Prange));
SumRate_RZF = zeros(1,length(Prange));
for num = 1:number
    num
    H = (randn(K,M)+1i*randn(K,M))*sqrt(0.5);
    for i = 1:length(Prange)
        P = Prange(i);
        W_ZF = ZF(H,P);
        W_MMSE = MMSE(H,P);
        W_RZF = RZF(H,P,path_loss);
        SumRate_ZF(i) = SumRate_ZF(i)+CalSumRate(H,W_ZF,path_loss);
        SumRate_MMSE(i) = SumRate_MMSE(i)+CalSumRate(H,W_MMSE,path_loss);
        SumRate_RZF(i) = SumRate_RZF(i)+CalSumRate(H,W_RZF,path_loss);
    end
end
SumRate_ZF = SumRate_ZF/number;
SumRate_MMSE = SumRate_MMSE/number;
SumRate_RZF = SumRate_RZF/number;