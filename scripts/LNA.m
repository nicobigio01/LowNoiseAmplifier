%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Nicolas Bigiotti - Alberto Moretti
%
%   LNA Design
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

% OPAMP characteristics
A_0_db = 40; %dB
A_0 = db2mag(A_0_db);
UGBW = 75e6; %Hz
PM =  80; %degree


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   OPAMP Transfer function calculation for cadence simulation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate the first dominant pole of the opamp TF
tau_PD = A_0 / (2*pi*UGBW);

tau_PS = 1 / (2*pi*((PM - 45) * 0.2 *UGBW + UGBW));

s = tf('s');
A = A_0 / ((1+s*tau_PD)*(1+s*tau_PS))

%controlSystemDesigner('bode',A);
bode(A)
f_highpass = 0.5e6;
LNA_gain_db =  10;
LNA_gain = db2mag(LNA_gain_db);

R_F = 100e3;
C_F = 1/(f_highpass*R_F)
C_S = C_F*db2mag(10)



T = A_0*(1+s*R_F*C_F)/(1+s*(C_F+C_S)*(R_F)*(1+s*tau_PD))
controlSystemDesigner('bode',T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   DCOP
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
