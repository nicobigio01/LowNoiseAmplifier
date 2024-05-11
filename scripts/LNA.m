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

% Useful options
bodeopts_mag = bodeoptions;

bodeopts_mag.FreqUnits = 'Hz';
bodeopts_mag.PhaseVisible = 'off';
bodeopts_mag.Grid = 'on';


bodeopts_phase = bodeoptions;

bodeopts_phase.FreqUnits = 'Hz';
bodeopts_phase.MagVisible = 'off';
bodeopts_phase.Grid = 'on';


addpath('/home/nicobigio01/cadence_vm/fromCadence/LabEle2');

currentImgPath = "";

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Design Constants
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Costants
k = 13.8e-24; % J/K
T = 300; % K


% OPAMP characteristics
A_0_db = 40; %dB
A_0 = db2mag(A_0_db);
UGBW = 75e6; %Hz
PM =  80; %degree

% LNA characteristics
f_hp_LNA = 0.5e6; %Hz
LNA_gain_db =  10; %dB
LNA_gain = db2mag(LNA_gain_db);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   OPAMP Transfer function evaluation for cadence simulation
%a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Calculate the first dominant pole of the opamp TF in order to respect
% UGBW
tau_A = A_0 / (2*pi*UGBW);

% Calculate the second pole from phase margin costraint
tau_B = tan(-deg2rad(-180+PM)-atan(tau_A*2*pi*UGBW))/(2*pi*UGBW);

% Adding a parametric right hand zero
tau_C = logspace(8,10,3);

% OPAMP Transfer Function
s = tf('s');
A = A_0 / ((1+s*tau_A)*(1+s*tau_B))

%controlSystemDesigner('bode',A);

[A_mag , A_phase , A_w] = bode(A);

h1 = figure();
semilogx(A_w./(2*pi) , mag2db(reshape( A_mag(1,1,:) , [1 96])));
hold on
xline(UGBW , "LineWidth" , 1.5 , "LineStyle" , "--" , "Color" , "red")
yline(0 , "LineWidth" , 1.5 , "LineStyle" , "--" , "Color" , "red")
text(UGBW - 50e6, -70 ,  "UGBW: 75 MHz" , "BackgroundColor" , "white" , "EdgeColor" , "red");
ylabel("Mag [dB]");
xlabel("Freq log[Hz]");
title("Modulo Modello OPAMP");
ylim([-80 45])
grid on;

h2 = figure();
semilogx(A_w./(2*pi) , reshape( A_phase(1,1,:) , [1 96]));
hold on;
xline(UGBW , "LineWidth" , 1.5 , "LineStyle" , "--" , "Color" , "red");
yline(-180+83.5 , "LineWidth" , 1.5 , "LineStyle" , "--" , "Color" , "red");
text(UGBW - 50e6, -150 ,  "PM: 83.5°" , "BackgroundColor" , "white" , "EdgeColor" , "red");
ylabel("Phase [°]");
xlabel("Freq [Hz]");
title("Fase Modello OPAMP");
grid on;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   DCOP
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Nothing to calculate!
% VCM si applied at node inp, for virtual ground principle inm is biased at
% same potential and so the output since the LNA configuration in DC is
% basically a buffer.




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Small Signal
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RF = logspace(3,5,3);

% By setting the R_F value and knowing that the pole frequency sets the
% highpass frequency corner we have:
CF = 1./(2*pi*f_hp_LNA.*RF);

% From the small signal gain at center frequency which is C_S/C_F we
% obtain:
CS = CF.*LNA_gain;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Noise
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_space = 2*pi*logspace(1,10,200);

% First we set the opamp noise to 5nV/sqrt(Hz)
v_nopamp = [5e-9 8e-9 10e-9]; %V/sqrt(Hz)



syms f R_F C_F C_S VNOPAMP

v_nr(R_F) = sqrt(4*k*T*R_F);

A_LNA(f , R_F , C_F , C_S , VNOPAMP) = ((2*pi*f*R_F*C_S)^2)/(1+(2*pi*f*R_F*C_F)^2);

VN_R(f , R_F , C_F , C_S , VNOPAMP) = (v_nr^2)/(1+(2*pi*f*R_F*C_F)^2);
VN_OPAMP(f , R_F , C_F , C_S , VNOPAMP) =(VNOPAMP^2*((2*pi*f*(C_S+C_F)*R_F)^2+1))/(1+(2*pi*f*R_F*C_F)^2);
VNIN_OPAMP(f , R_F , C_F , C_S , VNOPAMP) = VN_OPAMP / A_LNA;
VNIN_R(f , R_F , C_F , C_S , VNOPAMP) = VN_R / A_LNA;

VN(f , R_F , C_F , C_S , VNOPAMP) = VN_R+VN_OPAMP;
VNIN(f , R_F , C_F , C_S , VNOPAMP) = VN*(1+(2*pi*f*R_F*C_F)^2)/((2*pi*f*R_F*C_S)^2);


noise = [sqrt(double(VNIN(f_space,RF(1),CF(1),CS(1),v_nopamp(1)))) ; sqrt(double(VNIN(f_space,RF(2),CF(2),CS(2),v_nopamp(1)))) ; sqrt(double(VNIN(f_space,RF(3),CF(3),CS(3),v_nopamp(1))))];

h3 = figure();

loglog(f_space , noise(1,:));
grid on;
hold on;
loglog(f_space , noise(2,:));
loglog(f_space , noise(3,:));
title("R_F Contribution to input referred noise")
xlabel("f log[Hz]")
ylabel("v_{nin} log[V/sqrt(Hz)]")
legend("R_F=1kOhm" , "R_F=10kOhm" , "R_F=100kOhm")

noise = [sqrt(double(VNIN(f_space,RF(3),CF(3),CS(3),v_nopamp(1)))) ; sqrt(double(VNIN(f_space,RF(3),CF(3),CS(3),v_nopamp(2)))) ; sqrt(double(VNIN(f_space,RF(3),CF(3),CS(3),v_nopamp(3))))];


h4 = figure();

loglog(f_space , noise(1,:));
grid on;
hold on;
loglog(f_space , noise(2,:));
loglog(f_space , noise(3,:));
title("Miller OTA Contribution to input referred noise")
xlabel("f log[Hz]")
ylabel("v_{nin} log[V/sqrt(Hz)]")
legend("v_{n,OTA}=5nV/sqrt(Hz)" , "v_{n,OTA}=8nV/sqrt(Hz)" , "v_{n,OTA}=10nV/sqrt(Hz)")







GLoop = -A*(1+s*RF(3)*CF(3))/(1+s*RF(3)*(CF(3)+CS(3)));
[GLoop_mag , GLoop_phase , GLoop_w] = bode(GLoop);

h5 = figure();
semilogx(GLoop_w./(2*pi) , mag2db(reshape( GLoop_mag(1,1,:) , [1 110])));
hold on
xline(18.9E6 , "LineWidth" , 1.5 , "LineStyle" , "--" , "Color" , "red")
yline(0 , "LineWidth" , 1.5 , "LineStyle" , "--" , "Color" , "red")
text(18.9E6, -70 ,  "f: 18.9 MHz" , "BackgroundColor" , "white" , "EdgeColor" , "red");
ylabel("Mag [dB]");
xlabel("Freq log[Hz]");
title("G_{LOOP} - R_F=100kOhm, C_S=3pF, C_F=11pF");
ylim([-80 45])
grid on;

h6 = figure();
semilogx(GLoop_w./(2*pi) , reshape( GLoop_phase(1,1,:) , [1 110]));
hold on;
xline(18.9E6 , "LineWidth" , 1.5 , "LineStyle" , "--" , "Color" , "red");
yline(88.91 , "LineWidth" , 1.5 , "LineStyle" , "--" , "Color" , "red");
text(UGBW - 50e6, -150 ,  "PM: 88.4°" , "BackgroundColor" , "white" , "EdgeColor" , "red");
ylabel("Phase [°]");
xlabel("Freq [Hz]");
title("Fase Modello OPAMP");
grid on;


GLoop = -A*(1-s*6.28e-9)*(1+s*RF(3)*CF(3))/(1+s*RF(3)*(CF(3)+CS(3)));
[GLoop_mag , GLoop_phase , GLoop_w] = bode(GLoop);

h7 = figure();
semilogx(GLoop_w./(2*pi) , mag2db(reshape( GLoop_mag(1,1,:) , [1 110])));
hold on
xline(18.9E6 , "LineWidth" , 1.5 , "LineStyle" , "--" , "Color" , "red")
yline(0 , "LineWidth" , 1.5 , "LineStyle" , "--" , "Color" , "red")
text(18.9E6, -70 ,  "f: 18.9 MHz" , "BackgroundColor" , "white" , "EdgeColor" , "red");
ylabel("Mag [dB]");
xlabel("Freq log[Hz]");
title("G_{LOOP} - R_F=100kOhm, C_S=3pF, C_F=11pF");
ylim([-80 45])
grid on;

h8 = figure();
semilogx(GLoop_w./(2*pi) , reshape( GLoop_phase(1,1,:) , [1 110]));
hold on;
xline(18.9E6 , "LineWidth" , 1.5 , "LineStyle" , "--" , "Color" , "red");
yline(51.55 , "LineWidth" , 1.5 , "LineStyle" , "--" , "Color" , "red");
text(UGBW - 50e6, -150 ,  "PM: 88.4°" , "BackgroundColor" , "white" , "EdgeColor" , "red");
ylabel("Phase [°]");
xlabel("Freq [Hz]");
title("Fase Modello OPAMP");
grid on;

%controlSystemDesigner('bode',GLoop);


%T = A_0*(1+s*R_F*C_F)/(1+s*(C_F+C_S)*(R_F)*(1+s*tau_A))
%%controlSystemDesigner('bode',T);
%
%tmp = importdata("LabEle2_SystemDesignTest_FrequencyResponse.matlab");
%cadence_FreqRes_f = tmp.data(:,1);
%cadence_FreqRes_A = tmp.data(:,2);
%tmp = importdata("LabEle2_SystemDesignTest_SineOut.matlab");
%cadence_SineOut_t = tmp.data(:,1);
%cadence_SineOut_V = tmp.data(:,2);
%tmp = importdata("LabEle2_SystemDesignTest_SineIn.matlab");
%cadence_SineIn_t = tmp.data(:,1);
%cadence_SineIn_V = tmp.data(:,2);
%
%
%% Find 3dB frequency
%[Max_gain,index1]=max(cadence_FreqRes_A);
%freq_max=cadence_FreqRes_f(index1);
%Gain_3db=Max_gain-3;
%idxr = find(diff(sign(cadence_FreqRes_A-Gain_3db)))+(-1:1);
%
%for k = 1:size(idxr,1)
%    f_3db(k) = interp1(cadence_FreqRes_A(idxr(k,:)),cadence_FreqRes_f(idxr(k,:)),Gain_3db,'spline');
%end
%
%
%h3 = figure();
%
%semilogx(cadence_FreqRes_f , cadence_FreqRes_A)
%hold on;
%yline(max(cadence_FreqRes_A)-3 , "LineWidth" , 1.5 , "LineStyle" , "--" , "Color" , "red");
%xline(f_3db(1) , "LineWidth" , 1.5 , "LineStyle" , "--" , "Color" , "red");
%xline(f_3db(2) , "LineWidth" , 1.5 , "LineStyle" , "--" , "Color" , "red");
%text(f_3db(1) - 300e3, -30 ,  sprintf("f high pass: %.1f kHz" , f_3db(1)/1e3) , "BackgroundColor" , "white" , "EdgeColor" , "red");
%text(f_3db(2) -10e6, -30 ,  sprintf("f low pass: %.1f MHz" , f_3db(2)/1e6) , "BackgroundColor" , "white" , "EdgeColor" , "red");
%grid on;
%legend("LNA Frequency Response");
%title("LNA Frequency Response");
%ylabel("Mag [dB]");
%xlabel("Freq [Hz]");
%ylim([-35 15]);
%
%h4 = figure();
%plot(cadence_SineOut_t,cadence_SineOut_V);
%hold on;
%plot(cadence_SineIn_t,cadence_SineIn_V);
%grid on;
%legend("v_{out}" , "v_{in}");
%title("LNA Transient Analysis");
%ylabel("Voltage [V]");
%xlabel("Time [s]");
%xlim([0 10e-6]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Results summary
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('*** OPAMP Model Summary ***\n');
fprintf('tau_A: %f\n' , tau_A);
fprintf('tau_A: %f\n' , tau_B);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Saving and cleanup
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% saveas(h1,"../latex/design_review1_26032024/img/LabEle2_LNA_opamp_model_tf_mag.png");
% saveas(h2,"../latex/design_review1_26032024/img/LabEle2_LNA_opamp_model_tf_phase.png");
saveas(h3,"../latex/design_review1_26032024/img/LabEle2_LNA_opamp_model_VNIN_R_contr.png");
saveas(h4,"../latex/design_review1_26032024/img/LabEle2_LNA_opamp_model_VNIN_OTA_contr.png");