function [pulse_playedout_V,PulseProps]=calculatePulseProperties_HanSinc(amp_pulse_V,tau_pulse_s,noSidelobes)

%% define basic constants
% nPulseSteps = 1000;
% dt = tau_pulse_s/nPulseSteps; 
% t=linspace(-tau_pulse_s/2,tau_pulse_s/2,nPulseSteps);
%% define basic constants
Z0=50;%Ohm=V/A; impedance
dt = 100e-9;%100ns raster 
nPulseSteps = tau_pulse_s/dt;
t = linspace(-tau_pulse_s/2,tau_pulse_s/2,nPulseSteps);
%% define pulse shape
%BWT=1.6; %[1]
%Emp=0.8;% this parameter is set in the source code. it is meant to correct for imperfections
%noSidelobes=2
han_filter = (1+cos(2*pi*t/tau_pulse_s))/2;
%correctionFactor=2/pi/Emp/BWT*(noSidelobes+1); 
%waveform_sinc=sinc(pi*BWT*Emp*t/tau_pulse_s* correctionFactor);
waveform_sinc = sinc(2*t/tau_pulse_s*(noSidelobes+1));

pulse_playedout_V = amp_pulse_V * han_filter .* waveform_sinc;%[V]
PulseProps.pulseIntegral_Vs = sum(pulse_playedout_V.*dt);%unit [V*s]
PulseProps.pulseSqrIntegral_sqrVs = sum(pulse_playedout_V.^2 .*dt);%unit [V^2*s]

PulseProps.energy_VAs=PulseProps.pulseSqrIntegral_sqrVs /Z0;
PulseProps.meanPower_W=PulseProps.energy_VAs /tau_pulse_s;
