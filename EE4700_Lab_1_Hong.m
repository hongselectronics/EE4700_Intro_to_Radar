%Name: Jeremy Hong
%Date: 9/23/2015
%Class: EE 4700 Introduction to Radar
%Lab 1
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
P_R = 125000;          %Peak transmit power P_R = 125 kW
G_RT = 36;             %Mainbeam gain G_RT = 36 dBi
f_c = 6000000000;      %Carrier Frequency f_c = 6 GHz
G_sp = 30;             %Signal processing gain G_sp = 30       
F_R = 7;               %Receiver noise figure F_R = 7 dB
B_R = 30000000;        %Receiver processing bandwidth B_R = 30 MHz
L_R = 9;               %Total radar related losses L_R = 9 dB
half_power = -3;       %Antenna half-power half_power = -3 dB
B_theta = 1.5;         %Beamwidth B_theta = 1.5 degrees
scan = 36;             %Antenna scan rate scan = 36 degrees/second
P_d = 0.50;            %Probability of detection P_d = 0.50
P_fa = 10^-6;          %Probability of fals alarm P_fa = 10^-6
PRF = 1500;            %Pulse repetition frequency PRF = 1500 Hz
%% Target Characteristics
RCS = 0;               %Radar cross section RCS = 0 dBsm
%
%% Wavelength calculation
Wavelength = c/f_c;    %Wavelength meters (m)
%% Conversions to absolute values
G_RT = 10^(G_RT/10);
F_R = 10^(F_R/10);
L_R = 10^(L_R/10);
half_power = 10^(half_power/10);
RCS = 10^(RCS/10);

%% Calculations
%
%% The received single-pulse target signal power (S) at RRT_test in both Watts and dBW  
S =(P_R*(G_RT^2)*(Wavelength^2)*RCS)/(((4*pi)^3)*(R_RT_test^4)*L_R);
S_dB = 10*log10(S);
fprintf('\nThe received single-pulse target signal power is: %d Watts\n',S)
fprintf('The received single-pulse target signal power is: %3.4f dBW\n\n',S_dB)
%% The radar receiver thermal noise power (N) in both Watts and dBW  
N = F_R*k*T_0*B_R;
N_dB = 10*log10(N);
fprintf('\nThe radar receiver thermal noise power is: %d Watts\n',N)
fprintf('The radar receiver thermal noise power is: %3.4f dBW\n\n',N_dB)
%% Plot the received single-pulse target signal power and receiver thermal noise power vs. the radar-to-target range: 0 ~ 100 km (similar to Figure 2-6)
%Plot (1)
range1 = 0:1:100000;
%Signal Power
S_Plot = (P_R.*(G_RT).^2*Wavelength.^2.*RCS)./((4.*pi).^3.*(range1).^4.*L_R);
%Conversion to dB
S_Plot_dB = 10*log10(S_Plot);
plot(range1,S_Plot_dB)
hold on
%Noise Power 
N_Plot = F_R*k*T_0*B_R;
%Conversion to dB
N_dB_Plot = ones(1,100001)*10*log10(N_Plot);
plot(range1,N_dB_Plot)
%Legend
legend('Signal (S)','Thermal Noise (N)')
%Set Axis like plot 2-7
axis([0 100000 -150 -60])
%Labels
title('Single Pulse recieved target signal power vs Radar-to-Target range') 
ylabel('Power (dBW)')
xlabel('Radar-to-target range (m)')
%% The single-pulse target signal-to-noise ratio (S/N) at RRT_test in both absolute and dB 
SNR = S*G_sp/N;
SNR_dB = 10*log10(SNR);
fprintf('Single Pulse target signal-to-noise Ratio: %3.4f \n',SNR);
fprintf('Single Pulse target signal-to-noise Ratio: %3.4f dB \n\n',SNR_dB);
%% The signal-to-noise ratio required for detection (SNR_dt), detection threshold, for a Swerling Case 1 target in both absolute and dB 
SW_1_SNR =  (log(P_fa)/log(P_d)) - 1;
SW_1_SNR_dB = 10*log10(SW_1_SNR);
fprintf('SNR required for detection: %3.4f \n',SW_1_SNR);
fprintf('SNR required for detection: %3.4f dB \n\n',SW_1_SNR_dB);
%% The integration gain (GI) for a Swerling Case 1 target in both absolute and dB 
T_ill = B_theta/scan;
T_I = T_ill;
n_p = T_I*PRF;
G_I = n_p;
G_I_dB = 10*log10(G_I);
fprintf('Integration gain for a Swerling Case 1 target is: %3.4f \n',G_I);
fprintf('Integration gain for a Swerling Case 1 target is: %3.4f dB \n\n',G_I_dB);
%% The multiple pulse target signal‐to‐noise ratio (S/N)n at RRT_test in both absolute and dB 
SNR_n = G_I*SNR;
SNR_n_dB = 10*log10(SNR_n);
fprintf('The multiple pulse target signal-to-noise ratio is: %3.4f \n',SNR_n);
fprintf('The multiple pulse target signal-to-noise ratio is: %3.4f dB \n\n',SNR_n_dB);
%% Plot the multiple pulse target signal-to-noise ratio, radar detection threshold (SNRdt), and single-pulse target signal-to-noise ratio vs. the radar‐to‐target range: 0 ~ 100 km (similar to Figure 3-13)
figure(2)
range2 = 0:1:100000;
%Plot the SNR
SNR_Plot = (P_R.*G_RT.^2.*Wavelength.^2.*RCS*G_sp)./((4.*pi).^3.*range2.^4.*F_R.*k.*T_0.*B_R.*L_R);
SNR_Plot_dB = 10*log10(SNR_Plot);
plot(range2,SNR_Plot_dB)
hold on
%Plot the SNR_dt
SNR_dt = SW_1_SNR;
SNR_dt_Plot = ones(1,100001)*10*log10(SNR_dt);
plot(range2,SNR_dt_Plot,':')
%Plot SNR_n
SNR_n_Plot = (P_R.*G_RT.^2.*Wavelength.^2.*RCS*G_sp*G_I)./((4.*pi).^3.*range2.^4.*F_R.*k.*T_0.*B_R.*L_R);
SNR_n_Plot_dB = 10*log10(SNR_n_Plot);
plot(range2,SNR_n_Plot_dB,'--')
%Legend
legend('(S/N)','SNR_d_t','(S/N)_n')
%Set Axis like plot 3-13
axis([0 100000 -10 80])
%Labels
title('Multiple Pulse SNR') 
ylabel('SNR (dB)')
xlabel('Radar-to-target range (m)') 
%%The radar detection range (Rdt) in meters 
R_dt = ((P_R*G_RT^2*Wavelength^2*RCS*G_sp)/((4*pi)^3*SNR_dt*F_R*k*T_0*B_R*L_R))^(1/4);
fprintf('Radar detection range: %3.3f meters \n',R_dt);

