%% Capacity Characterization for Intelligent Reflecting Surface Aided MIMO Communication
clc;clear;close all;
Nt = 64;Nr = 32;M = 128;
sigma_q = 1;
SNR_range = -10:5:20;
C_iter = [];
iter_max = 1000;
number = 100;
RISflag = 1;
C_all = zeros(length(SNR_range),number);
for snr = 1:length(SNR_range)
    SNR = SNR_range(snr);
    Pt = 10^(SNR/10);
    for num = 1:number
        num
        H = sqrt(1/2)*(randn(Nr,Nt)+1i*randn(Nr,Nt));
        T = sqrt(1/2)*(randn(M,Nt)+1i*randn(M,Nt));
        R = sqrt(1/2)*(randn(Nr,M)+1i*randn(Nr,M));
        alpha = exp(1i*2*pi*rand(M,1));
        for iter = 1:iter_max
            if RISflag == 0
                alpha = zeros(M,1);
            end
            H_tilde = H+R*diag(alpha)*T;
            [~,Lambda,V] = svd(H_tilde);
            D = rank(H_tilde);
            V_tilde = V(:,1:D);
            Lambda_tilde = Lambda(1:D,1:D);
            L = D;
            temp1 = sigma_q./diag(Lambda_tilde^2);
            while(1)
                level = (Pt+sum(temp1(1:L)))/L;%level = 1/p0
                if level>=temp1(L)
                    break
                end
                L = L - 1;
            end
            p = [level - temp1(1:L);zeros(D-L,1)];
            Q = V_tilde*diag(p)*V_tilde';
            C = sum(log2(1+diag(Lambda_tilde^2).*p/sigma_q));
            if RISflag == 0
                break;
            end
%             C1 = log2(real(det(eye(Nr)+1/sigma_q*H_tilde*Q*H_tilde')));
            C_iter = [C_iter,C];
            if length(C_iter)>10 && abs((C_iter(end)-C_iter(end-1))/C_iter(end))<1e-3
                break
            end
            [U_Q,SIGMA_Q] = eig(Q);
            H_prime = H*U_Q*SIGMA_Q^(1/2);
            T_prime = T*U_Q*SIGMA_Q^(1/2);
            T_prime = T_prime';
            for m = 1:M
                rm = R(:,m);
                tm_prime = T_prime(:,m);
                temp2 = zeros(Nr,Nt);
                for i = [1:m-1,m+1:M]
                    ri = R(:,i);
                    ti_prime = T_prime(:,i);
                    temp2 = temp2 + alpha(i)*ri*ti_prime';
                end
                Am = eye(Nr)+1/sigma_q*(H_prime+temp2)*(H_prime+temp2)'+1/sigma_q*(rm*tm_prime')*(rm*tm_prime')';
                Bm = 1/sigma_q*rm*tm_prime'*(H_prime+temp2)';
                C3 = log2(real(det(Am+(alpha(m)*Bm)+(alpha(m)*Bm)')));
                temp3 = trace(Am\Bm);
                if temp3 ~= 0 
                    alpha(m) = exp(-1i*angle(temp3));
                else
                    alpha(m) = 1;
                end
            end
        end
        C_all(snr,num) = C;
    end
end
% save('RISMIMO.mat');
save('NonRISMIMO.mat');
load('RISMIMO.mat');
plot(SNR_range,mean(C_all,2),'-k','linewidth',1.5);hold on;
load('NonRISMIMO.mat');
plot(SNR_range,mean(C_all,2),'-r','linewidth',1.5);hold on;
xlabel('SNR(dB)');
ylabel('rate(bps/Hz)');
legend('with RIS','without RIS');
grid on;