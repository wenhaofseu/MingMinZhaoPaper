clc;clear;close all;
rng(0);
M = 6;% number of AP antennas
N = 40;% number of IRS elements
K = 4;% number of users
sigma_q_k = 10^(-80/10-3); % noise power:-80dBm
%% Rician factor
betaIU = 10^(5/10); %IRS-User
betaAI = 10^(5/10); %AP-IRS
betaAU = 10^(-5/10); %AP-User
%% location info
dv = 2;% AP loaction parameter(m)
d0 = 50; % IRS location parameter(m)
pAP = [dv,0,0]; %AP location
pRIS = [0,d0,3]; %IRS location
pUser = [[-3,50-3]+6*rand(3,2),[0;0;0]]; %User location
%% correlation coefficient
r_r = 0.8; %PHI_r
r_d = 0.8; %PHId
r_rk = 0.8; %PHI_rk
%% correlation matrix
PHI_r = expCorModel(r_r,N);
PHId = expCorModel(r_d,M);
PHI_rk = expCorModel(r_rk,N);
%% deterministic component of channel


