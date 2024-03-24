%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Nicolas Bigiotti - Alberto Moretti
%
%   Miller OpAmp design
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Cleaning workspace
clc
clear
close all

%addpath('/home/nicobigio01/cadence_vm/fromCadence/lab_ele1_cadence_ex2');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Design Constants
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

q = 1.6e-19; % C
k = 13.8e-24; % J/K
T = 300; % K
epsilon_0 = 8.89e-12; % F/m
epsilon_si = 11.3 * epsilon_0;
epsilon_ox = 3.9 * epsilon_0;
N_A = 4e23; % m^(-3)
t_ox = 4e-9; % m
mu_n = 0.0348; % m^2/(V*s)
V_builtin = 0.6; % V
V_TH = 0.55; % V
VDD = 1.8; %V
lambda_prime = 5e-8;


K_n = 0.5 * mu_n * epsilon_ox / t_ox;
%K_p = 35.6e-6;
K_p = K_n / 3;
C_ox = epsilon_ox / t_ox;
C_L = -1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Design Parameter
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DCOP
V_OV1 = -1; V_OV2 = -1; V_OV3 = -1; V_OV4 = -1; V_OV6 = -1;
I_1 = -1; I_2 = -1; I_3 = -1; I_4 = -1; I_6 = -1;
gm1 = -1; gm2 = -1; gm3 = -1; gm4 = -1; gm6 = -1;
rds1 = -1; rds2 = -1; rds3 = -1; rds4 = -1; rds6 = -1;

% MOS Sizing
W1 = -1; W2 = -1; W3 = -1; W4 = -1; W6 = -1; 
L1 = -1; L2 = -1; L3 = -1; L4 = -1; L6 = -1; 

% Compensation 
C_C = -1; R_C = -1;

% Total Bias Current
I_A = -1; I_B = -1;




A_v_db = 40; %dB
A_v = db2mag(A_v_db);
UGBW = 75e6; %Hz
PM = 80; %degree
PSD = 10e-9; %V/sqrt(HZ)

%SR=10e7; %V/s

A_1 = 10;
A_2 = 10;


V_OV1 = 0.4;

gm1 = (32*k*T)/(3*PSD*PSD)

I_1 = 0.5 * V_OV1 * gm1

I_2 = I_1;
I_3 = I_1;
I_4 = I_1;

I_A = 2 * I_1;

L1 = 2 * lambda_prime * I_1 * A_1 / gm1



% L_12 = A_1 *lambda_prime * 0.2
% W_12 = I_A*L_12 / (K_n*V_OV_12*V_OV_12)
% 
% L_34 = L_12;
% W_34 = K_n / K_p * W_12
% 
% 
% W_56 = 2*I_A*L_12 / (K_p*V_OV_12*V_OV_12)
% 
% 
% I_B = 0.5* I_A * W_56 * L_34 / (W_34*L_34)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Da pdf
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Saving and cleanup
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%saveas(f1,"../latex/img/lab_ele1_pract_ex3_static.png");
%saveas(f2,"../latex/img/lab_ele1_pract_ex3_mod_inv.png");
%saveas(f3,"../latex/img/lab_ele1_pract_ex3_phase_inv.png");
%saveas(f4,"../latex/img/lab_ele1_pract_ex3_mod_rlc.png");
%saveas(f5,"../latex/img/lab_ele1_pract_ex3_phase_rlc.png");
