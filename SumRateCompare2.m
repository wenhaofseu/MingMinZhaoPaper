clc;clear;close all;
number = 10000;%number of Monter Carlo
table = 2:1:10;
P_dB= 0;
P = 10^(P_dB/10);
SumRate_ZF = zeros(1,length(table));
SumRate_MMSE = zeros(1,length(table));
SumRate_RZF = zeros(1,length(table));
for i = 1:length(table)
    M = table(i);
    K = M;
    for num = 1:number
        path_loss = 110-130+unifrnd(-5,5,1,K);%noise power:-80dBm;average path loss:130dB
        rho = 10.^(path_loss/10);
        H = (randn(K,M)+1i*randn(K,M))*sqrt(0.5);
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