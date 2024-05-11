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
lambda_prime_n = 5e-8;  % Must fix these two values
lambda_prime_p = 5e-8;  
K_n = 0.5 * mu_n * epsilon_ox / t_ox;
K_p = K_n / 3;
C_ox = epsilon_ox / t_ox;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Design Specs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f_hi = 500e3; %Hz
f_lo = 25e6; %Hz
LNA_gain_db =  10; %dB
LNA_gain = db2mag(LNA_gain_db);

%%% Miller Opamp gain, phase margin and UGBW
A_0_db = 40; %dB
A_0 = db2mag(A_0_db);
UGBW = f_lo * LNA_gain; %Hz    % 75MHz from design specs, needs to be a little higher
                    % to accomodate the 25MHz low pass frequency of the full LNA
PM = 85; %deg   % 80 deg from specs

% Load capacitor and resistance
% These values are just assumptions
C_L = 1e-12; %F
R_L = 1e9; %Ohm

C_F = C_L;

R_F = 1./(2*pi*f_hi.*C_F);

% From the small signal gain at center frequency which is C_S/C_F we
% obtain:
C_S = C_F.*LNA_gain;

%Input referred noise spec
VIRN = 5e-9; %V/sqrt(Hz)

% Poles estimation from two dominant pole model of opamp
tau_A = A_0 / (2*pi*UGBW);
tau_B = tan(-deg2rad(-180+PM)-atan(tau_A*2*pi*UGBW))/(2*pi*UGBW);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Design Procedure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Defining dummy variables
% DCOP
V_OV1 = -1; V_OV2 = -1; V_OV3 = -1; V_OV4 = -1; V_OV5 = -1; V_OV6 = -1; V_OV7 = -1;
I1 = -1; I2 = -1; I3 = -1; I4 = -1; I5 = -1; I6 = -1; I7 = -1; I8 = -1;
gm1 = -1; gm2 = -1; gm3 = -1; gm4 = -1; gm5 = -1; gm6 = -1; gm7 = -1; gm8 = -1;
rds1 = -1; rds2 = -1; rds3 = -1; rds4 = -1; rds5 = -1; rds6 = -1; rds7 = -1; rds8 = -1;

% MOS Sizing
W1 = -1; W2 = -1; W3 = -1; W4 = -1; W5 = -1; W6 = -1; W7 = -1; W8 = -1;
L1 = 500e-9; L2 = 500e-9; L3 = 500e-9; L4 = 500e-9; L5 = 1500e-9; L6 = 500e-9; L7 = 1500e-9; L8 = 1500e-9;

% Compensation 
C_C = -1; R_C = -1;

% Total Bias Current
IA = -1; IB = -1;


%%% Design

IREF = 100e-6; % Fix a reasonable value for the reference current

% Keep C_C 2-4 times greater than C_L to comply to assumptions
C_C = 3*C_L;

% Evaluate gm1 constraint from bandwidth and noise requirements
gm1_noise = (16*k*T)/(3*VIRN^2); % gm1 needed for noise constraint 
gm1_bw = 2 * pi * UGBW * C_C; % gm1 needed for bandwidth constraint 

gm1 = gm1_bw*1.30; % Boost gm a little
gm2 = gm1;

% Fix overdrive of M1, lower means less current but also larger devices
V_OV1 = 0.15; %V
V_OV2 = V_OV1;

% Fix first stage currents
I1 = 0.5 * V_OV1 * gm1;
I2 = I1;
I3 = I1;
I4 = I1;
IA = 2 * I1;
I5 = IA;

%%% Second Stage

% Evaluate gm from phase margin costraint
gm6 = 1.3 * C_L / tau_B;

% Fix M6 overdrive remebering that is the same of M3 and M4
V_OV6 = 0.2;
V_OV3 = V_OV6;
V_OV4 = V_OV6;

% Calculate all currents
I6 = 0.5 * V_OV6 * gm6;
IB = I6;
I7 = IB;

% Fix all current mirror MOS overdrives
V_OV5 = 0.25;
V_OV7 = 0.25;
V_OV8 = 0.25;
I8 = IREF;


% Evaluate all Ws
W1 = (I1*L1) / ( K_n * V_OV1^2 );
W2 = (I2*L2) / ( K_n * V_OV2^2 );
W3 = (I3*L3) / ( K_p * V_OV3^2 );
W4 = (I4*L4) / ( K_p * V_OV4^2 );
W5 = (I5*L5) / ( K_n * V_OV5^2 );
W6 = (I6*L6) / ( K_p * V_OV6^2 );
W7 = (I7*L7) / ( K_n * V_OV7^2 );
W8 = (I8*L8) / ( K_n * V_OV8^2 );

% Evaluate all rds'
rds1 = (L1)/(lambda_prime_n * I1);
rds2 = (L2)/(lambda_prime_n * I2);
rds3 = (L3)/(lambda_prime_n * I3);
rds4 = (L4)/(lambda_prime_n * I4);
rds5 = (L5)/(lambda_prime_n * I5);
rds6 = (L6)/(lambda_prime_n * I6);
rds7 = (L7)/(lambda_prime_n * I7);
rds8 = (L8)/(lambda_prime_n * I8);


% Evaluate the input referred noise of the first stage without approximations
VIRN_tot =  sqrt((16*k*T)/(3*(gm1+gm6)));

% Evaluate the value of the zero nulling resistor
R_C=1/gm6;

%Evaluate gain
A_expected = (gm1*rds1*rds2*gm6*rds6*rds7)/((rds1+rds2)*(rds6+rds7));
A_expected_db = mag2db(A_expected);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Resume
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fprintf("Final Resume\nC_F: %.1f pF\nC_S: %.1f pF\nR_F: %.1f kOhm\nExpected Gain: %.1f dB\n\n" ,C_F*1e12 , C_S*1e12 , R_F/1e3);

fprintf("Final Resume\nVNIN_tot: %.1f nV/sqrt(Hz)\nC_C: %.1f pF\nR_C: %.1f Ohm\nExpected Gain: %.1f dB\n\n" , VIRN_tot*1e9 , C_C*1e12 , R_C , A_expected_db);


fprintf("M1 sizing resume:\nI1 = %.1f uA\nV_OV1 = %.2f V\nW1 = %.1f um\nL1 = %.1f nm\ngm1 = %.1f mA/V\nrds1 = %.1f kOhm\n\n\n" , ...
    I1*1e6 , ...
    V_OV1 , ...
    W1*1e6 ,...
    L1*1e9 , ...
    gm1*1e3 , ...
    rds1*1e-3 ...
    );

fprintf("M2 sizing resume:\nI2 = %.1f uA\nV_OV2 = %.2f V\nW2 = %.1f um\nL2 = %.1f nm\ngm2 = %.1f mA/V\nrds2 = %.1f kOhm\n\n\n" , ...
    I2*1e6 , ...
    V_OV2 , ...
    W2*1e6 ,...
    L2*1e9 , ...
    gm2*1e3 , ...
    rds2*1e-3 ...
    );

fprintf("M3 sizing resume:\nI3 = %.1f uA\nV_OV3 = %.2f V\nW3 = %.1f um\nL3 = %.1f nm\ngm3 = %.1f mA/V\nrds3 = %.1f kOhm\n\n\n" , ...
    I3*1e6 , ...
    V_OV3 , ...
    W3*1e6 ,...
    L3*1e9 , ...
    gm3*1e3 , ...
    rds3*1e-3 ...
    );


    

fprintf("M4 sizing resume:\nI4 = %.1f uA\nV_OV4 = %.2f V\nW4 = %.1f um\nL4 = %.1f nm\ngm2 = %.1f mA/V\nrds4 = %.1f kOhm\n\n\n" , ...
    I4*1e6 , ...
    V_OV4 , ...
    W4*1e6 ,...
    L4*1e9 , ...
    gm4*1e3 , ...
    rds4*1e-3 ...
    );

fprintf("M5 sizing resume:\nI5 = %.1f uA\nV_OV5 = %.2f V\nW5 = %.1f um\nL5 = %.1f nm\ngm5 = %.1f mA/V\nrds5 = %.1f kOhm\n\n\n" , ...
    I5*1e6 , ...
    V_OV5 , ...
    W5*1e6 ,...
    L5*1e9 , ...
    gm5*1e3 , ...
    rds5*1e-3 ...
    );

fprintf("M6 sizing resume:\nI6 = %.1f uA\nV_OV6 = %.2f V\nW6 = %.1f um\nL6 = %.1f nm\ngm6 = %.1f mA/V\nrds6 = %.1f kOhm\n\n\n" , ...
    I6*1e6 , ...
    V_OV6 , ...
    W6*1e6 ,...
    L6*1e9 , ...
    gm6*1e3 , ...
    rds6*1e-3 ...
    );

fprintf("M7 sizing resume:\nI7 = %.1f uA\nV_OV7 = %.2f V\nW7 = %.1f um\nL7 = %.1f nm\ngm7 = %.1f mA/V\nrds7 = %.1f kOhm\n\n\n" , ...
    I7*1e6 , ...
    V_OV7 , ...
    W7*1e6 ,...
    L7*1e9 , ...
    gm7*1e3 , ...
    rds7*1e-3 ...
    );

fprintf("M8 sizing resume:\nI8 = %.1f uA\nV_OV8 = %.2f V\nW8 = %.1f um\nL8 = %.1f nm\ngm8 = %.1f mA/V\nrds8 = %.1f kOhm\n\n\n" , ...
    I8*1e6 , ...
    V_OV8 , ...
    W8*1e6 ,...
    L8*1e9 , ...
    gm8*1e3 , ...
    rds8*1e-3 ...
    );

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
