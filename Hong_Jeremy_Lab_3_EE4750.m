%Name: Jeremy Hong
%Date: 10/28/2015
%Class: EE 4700 Introduction to Radar
%Lab 3
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
R_RT_test = 60000;     %Test radar-to-target range R_RT_test = 60 km
%% Radar Characteristics
del_t = 0.1*10^-3;     %Time delay TX/RX del_t = 0.1 mSec
f_cr = 9.000001*10^9;  %Recieved Carrier Frequency f_cr = 9.000001 GHz
f_c = 9e9;             %Carrier Frequency f_c = 9 GHz
tau = 0.5*10^-6;       %Transmitted pulse width tau = 5 uSec
PRF = 7500;            %Pulse repetition frequency PRF = 7500 Hz
T_I = 25*10^-3;        %Integration time (pulse burst duration) T_I = 25 mSec
N_phi = 100;           %Phase modulation pulse N_phi = 100
S_N = -3;              %Single pulse target signal-to-noise ration S/N = -3 dB
theta_3dB = 1.5;       %Circular antenna with half-power beamwidth theta_3dB = 1.5 degrees
ohm = 720;             %Total solid angle scan coverage omega = 720 degrees 
%% Wavelength calculation
Wavelength = c/f_cr;   %Wavelength meters (m)
%% Calculations
%% Conversions to absolute values and radians
S_N = 10^(S_N/10);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculations
%
%% The radar-to-target slant range (R) in meters
R = (c*del_t)/2;
fprintf('\n1.) The radar-to-target slant range is: %d meters\n',R)
%
%% The range resolution without pulse compression (del_R) in meters
del_R = (c*tau)/2;
fprintf('\n2.) The range resolution without pulse compression is: %d meters\n',del_R)
%
%% The range resolution with pulse compression (del_Rpc) in meters
del_R_pc = c/(2*(N_phi/tau));
fprintf('\n3.) The range resolution with pulse compression is: %3.4f meters\n',del_R_pc)
%
%% Plot the range resolution with pulse compression vs. pulse compression modulation BW (B_pc):
%% 0.1 to 100 MHz (similar to Figure 5-14)
fprintf('\n4.) See Figure 1\n');
B_pc = 100000:100:100000000;
del_R_pc_plot = c./(2.*B_pc);
figure()
loglog(B_pc,del_R_pc_plot);
grid on
title('Range Resolution as a function of Modulation Bandwidth')
xlabel('Modulation Bandwidth (MHz)')
ylabel('Range Resolution (m)')
%
%% The pulse compression signal processing gain (G_pc)
B_R = 1/tau;
B_pc = B_R;
G_pc = tau*B_pc;
fprintf('\n5.) The pulse compression signal processing gain is: %d absolute\n',G_pc);
%
%% The unambiguous range (R_u) in meters
R_u = c/(2*PRF);
fprintf('\n6.) The unambiguous range is: %d meters\n',R_u);
%
%% Plot the unmbiguous range vs. pulse repetition frequency (PRF): 100 to 10000 Hz (similar to figure 5-16)
fprintf('\n7.) See Figure 2\n');
PRF_plot = 100:1:10000;
R_u_plot = c./(2.*PRF_plot);
figure()
semilogx(PRF_plot,R_u_plot/1000);
grid on
title('Unambiguous Range a function of PRF')
xlabel('Pulse Repetition Frequency (PRF) (Hz)')
ylabel('Unambiguous Range (km)')
%
%% The dead zone (d_z) in meters
d_z = (c*tau)/2;
fprintf('\n8.) The dead zone is: %d meters\n',d_z);
%
%% The radar-to-target range rate (R_dot) in meters/second
f_d = f_cr - f_c; 
R_dot = (-f_d*Wavelength)/2;
fprintf('\n9.) The radar-to-target range rate is: %d meters/second\n',R_dot);
%
%% The range rate resolution (del_R_dot) in meters/second
del_R_dot = Wavelength/(2*T_I);
fprintf('\n10.) The range rate resolution is: %d meters/second\n',del_R_dot);
%
%% The unambiguous range rate (R_dotu) in meters/second
R_dotu = (PRF*Wavelength)/4;
fprintf('\n11.) The unambiguous range rate is: %d meters/second\n',R_dotu);
%
%% The angle resolution (del_theta) in degrees
del_theta = theta_3dB;
fprintf('\n12.) The angle resolution is: %3.1f degrees\n',del_theta);
%% The cross-range resolution (del_CR) in meters
del_CR = R*(theta_3dB*pi)/180;
fprintf('\n13.) The cross-range resolution is: %3.5f meters\n',del_CR);
%
%% The measurement accuracies:
sig_R = del_R/sqrt(2*G_pc*S_N);
sig_R_pc = del_R_pc/sqrt(2*G_pc*S_N);
sig_R_dot = del_R_dot/sqrt(2*G_pc*S_N);
sig_theta = del_theta/sqrt(2*G_pc*S_N);
fprintf('\n14.) The Measurement accuracies:\n');
fprintf('     Range without pulse compression: %3.4f meters\n',sig_R);
fprintf('     Range with pulse compression: %3.4f meters\n',sig_R_pc);
fprintf('     Range rate: %3.4f meters/second\n',sig_R_dot);
fprintf('     Angle: %3.4f degrees\n',sig_theta);
%
%%
fprintf('\n15.) Radar Resolution Cell \n');
n_rg = floor(R_u/del_R);
n_rg_pc = floor(R_u/del_R_pc);
del_f_d = 1/T_I;
n_df = floor(PRF/del_f_d);
n_b = ohm/theta_3dB;
fprintf('     Number of range gates without pulse compression: %d\n',n_rg);
fprintf('     Number of range gates with pulse compression: %d\n',n_rg_pc);
fprintf('     Number of Doppler filters: %d\n',n_df);
fprintf('     Number of antenna beam positions: %d \n',n_b);
