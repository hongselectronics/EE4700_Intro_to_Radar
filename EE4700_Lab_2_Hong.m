%Name: Jeremy Hong
%Date: 10/14/2015
%Class: EE 4700 Introduction to Radar
%Lab 2
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
D = 0.6;               %Antenna Diameter D = 0.6 meters
f_c = 6000000000;      %Carrier Frequency f_c = 6 GHz
e_space = 0.5;         %Element spacing relative to a wavelength d/wavelength = 0.5
rho = 0.60;            %Efficiency rho = 0.60
theta_0 = -30;         %Array antenna beam steering angle theta_0 = -30 degrees
%% Wavelength calculation
Wavelength = c/f_c;    %Wavelength meters (m)
%% Calculations
A = pi*(D/2)^2;        %Physical area of the antenna, square meters
%% The antenna mainbeam Gain (G) in both absolute and dBi
G = (4*pi*rho*A)/Wavelength^2;
G_dBi = 10*log10(G);
fprintf('The antenna mainbeam gain (G) is: %3.4f absolute, %3.4f dBi\n', G, G_dBi);
%
%% The half-power (-3 db) beamwidth (theta_3dB) in degrees
theta_3dB = 51*Wavelength/D;
fprintf('\nThe half-power (-3 dB) beamwidth is: %3.4f degrees\n',theta_3dB);
%
%% The null-to-null beamwidth (theta_nn) in degrees
theta_nn = 115*Wavelength/D;
fprintf('\nThe null-to-null beawidth is: %3.4f degrees\n',theta_3dB);
%
%% The antenna gain in the direction of the first sidelobe (G_1sl) in both absolute and dBi
% From table 4-1 the uniform illumination function has a first sidelobe
% level of -13.2 dBi
G_1sl_dBi = -13.2;
G_1sl = 10^(G_1sl_dBi/10);
fprintf('\nThe antenna gain in the direction of the first sidelobe is: %3.4f absolute, %3.4f dBi\n', G_1sl, G_1sl_dBi)
%
%% Plot the normalized antenna gain pattern vs. angle -90 to +90 degrees
theta = -pi/2:0.001:pi/2;
G_theta = (sin((pi*D*sin(theta))./Wavelength)./((pi*D*sin(theta))./Wavelength)).^2;
plot((180/pi)*theta, G_theta)
title('Normalized Antenna Gain Pattern vs. Angle')
xlabel('Angle (Degrees)')
ylabel('Normalized Gain')
%% The incremental phase shift between elements (delta_phi) necessary to steer the beam to an
%% angle (theta_0) in both radians and degrees.
delta_phi = 2*pi*(e_space)*sin(theta_0*pi/180);
delta_phi_deg = delta_phi*(180/pi);
fprintf('\nThe incremental phase shift between elements necessary to steer the beam to an angle is: %3.4f radians, %3.1f degrees\n',delta_phi,delta_phi_deg)
%
%%The array antenna mainbeam gain (G) when the beam is steered to an angle
%%(theta_0) in both absolute and dBi.
G_1 = G*cos(theta_0*pi/180);
G_1_dBi = 10*log10(G_1);
fprintf('\nThe array antenna mainbeam gain (G) when the beam is steered to angle (theta_0) is: %3.4f absolute, %3.4f dBi\n',G_1,G_1_dBi)
%
%% The half-power beamwidth (theta_3dB) when the beam is steered to an angle (theta_0) in degrees
theta_3dB1 = theta_3dB/cos(theta_0*pi/180);
theta_3dB1_deg = theta_3dB1*(180/pi);
fprintf('\nThe half-power beamwidth (theta_3dB) when the beam is steered to an angle is: %3.4f degrees\n',theta_3dB1_deg)
%
%% Plot the normalized array factor (G_a) vs. angle, -90 to +90 degrees, when the beam is steered to an angle (theta_0)


