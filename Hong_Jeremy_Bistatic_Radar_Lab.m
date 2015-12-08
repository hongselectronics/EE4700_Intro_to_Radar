%Name: Jeremy Hong
%Date: 12/08/2015
%Class: EE 4700 Introduction to Radar
%Bistatic Radar Software Lab Experiments
%
% Copyright notice & terms of use available at: 
% https://github.com/hongselectronics/EE4700_Intro_to_Radar/blob/master/LICENSE.md
% View code revision history here: 
%https://github.com/hongselectronics/EE4700_Intro_to_Radar
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Clean up
clear
clc
close all
format short
%
% 
%% Constants
c = 3*10^8;            %Speed of Light 3e8 m/sec
k = 1.38*10^-23;       %Boltzmann's constant 1.38e-23 watt-second/Kelvin
T_0 = 290;             %Standard reference temperature T_0 = 290 Kelvin
k_r = 4/3;             %refraction factor k_r = 4/3
R_E = 6371e3;          %Radius of Earth R_E = 6371 km
%% Bistatic Radar Characteristics
P_R = 50e3;            %Peak transmit power P_R = 50 kW
G_Rt_T = 0;            %Tx Ant Gain in the direction of the target G_Rt_t = 0 dBi
G_Rr_T = 10;           %Rx Ant Gain in the direction of the target G_Rr_t = 10 dBi
f_c = 100e6;           %Carrier Frequency f_c = 100 MHz
T_sp = 0.25;           %Signal processing interval T_sp = 0.25 s      
F_R = 6;               %Receiver noise figure F_R = 6 dB
B_R = 200e3;           %Receiver processing bandwidth B_R = 200 kHz
L_R = 7;               %Total radar related losses L_R = 7 dB
h_Rt = 100;            %Tx height h_Rt = 100 m
h_Rr = 10;             %Rx height h_Rr = 10 m
%% Target Characteristics
sig_b = 5;             %Bistatic RCS sigma_b = 5 m^2
h_T = 2000;            %Height h_T = 2000 m
R1 = 200e3;            %Tx-to-Target range R1 = 200 km
R2 = 100e3;            %Rx-to-Target range R2 = 100 km
%
%% Wavelength calculation
Wavelength = c/f_c;    %Wavelength meters (m)
%% Conversions from dB to Absolute values
G_Rt_T = 10^(G_Rt_T/10);
G_Rr_T = 10^(G_Rr_T/10);
F_R = 10^(F_R/10);
L_R = 10^(L_R/10);
%
%% Calculations
%
%% 1. The signal processing gain (G_sp)
G_sp = T_sp * B_R;
fprintf('1. The signal processing gain is: %d \n',G_sp);
%
%% 2. The processed target signal-to-noise ratio (S/N)_p in both absolute and decibels
S = P_R * G_Rt_T * G_Rr_T * Wavelength^2 * sig_b * G_sp;
N = (4*pi)^3 * R1^2 * R2^2 * F_R * k * T_0 * B_R * L_R;
SNR_p = S/N;
SNR_p_dB = 10*log10(SNR_p);
fprintf('2. The processed target SNR is: %3.4f absolute\n',SNR_p)
fprintf('   The processed target SNR is: %3.4f dB\n',SNR_p_dB)
%
%% 3. The bistatic detection range product (R1R2) in both square meters 
%% (m^2) and square kilometers (km^2)
R1R2_m = R1*R2;
R1R2_km = R1R2_m/1000;
fprintf('3. The bistatic detection range product is: %3.1E m^2 or %3.1E km^2\n',R1R2_m,R1R2_km);
%% 4. Plot the bistatic detection range coverage as a function of down-range and cross-range.
y = -300e3:0.1:300e3;
x = -300e3:0.1:300e3;
R1R2_sq = [(x-200e3./2).^2 + y.^2].*[(x+200e3./2).^2 + y.^2];
y1 = .5.*sqrt(-4*x.^2-200e3^2+4*sqrt(x.^2*200e3^2+R1R2_sq));
plot(R1R2_sq,y1,-R1R2_sq,y1,R1R2_sq,-y1,-R1R2_sq,-y1)
hold on
R1R2_sq = [(x-300e3./2).^2 + y.^2].*[(x+300e3./2).^2 + y.^2];
y1 = .5.*sqrt(-4*x.^2-300e3^2+4*sqrt(x.^2*300e3^2+R1R2_sq));
plot(R1R2_sq,y1,-R1R2_sq,y1,R1R2_sq,-y1,-R1R2_sq,-y1)
R1R2_sq = [(x-400e3./2).^2 + y.^2].*[(x+400e3./2).^2 + y.^2];
y1 = .5.*sqrt(-4*x.^2-400e3^2+4*sqrt(x.^2*400e3^2+R1R2_sq));
plot(R1R2_sq,y1,-R1R2_sq,y1,R1R2_sq,-y1,-R1R2_sq,-y1)

