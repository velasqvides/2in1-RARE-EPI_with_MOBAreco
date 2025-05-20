function RFpulseProps=calculatePulseProperties_HanSinc_BWT(amp_pulse_V,tau_pulse_s,BWT)

%% define basic constants
% nPulseSteps = 1000;
% dt = tau_pulse_s/nPulseSteps; 
% t=linspace(-tau_pulse_s/2,tau_pulse_s/2,nPulseSteps);
%% define basic constants
Z0=50;%Ohm=V/A; impedance
dt = 1e-9;%1ns raster 
nPulseSteps = tau_pulse_s/dt;
t = linspace(-tau_pulse_s/2,tau_pulse_s/2,nPulseSteps);
%% define pulse shape
%Emp=0.8;% this parameter is set in the source code. it is meant to correct for imperfections
han_filter = (1+cos(2*pi*t/tau_pulse_s))/2;
%correctionFactor=2/pi/Emp/BWT*(noSidelobes+1); 
waveform_sinc=sinc(BWT*t/tau_pulse_s);
pulse_playedout_V = amp_pulse_V * han_filter .* waveform_sinc;%[V]

RFpulseProps.sourceOfRFpulse='calculated hanning-filtered sinc';
RFpulseProps.transmitBandwidthTimeProd=BWT;

RFpulseProps.peakAmplitude_V=amp_pulse_V;
RFpulseProps.pulse_duration_s=tau_pulse_s;

RFpulseProps.pulse_playedout_V = pulse_playedout_V ;
RFpulseProps.rastertime_s=dt;

RFpulseProps.pulseIntegral_Vs = sum(pulse_playedout_V.*dt);%unit [V*s]
RFpulseProps.pulseSqrIntegral_sqrVs = sum(pulse_playedout_V.^2 .*dt);%unit [V^2*s]
RFpulseProps.energy_VAs=RFpulseProps.pulseSqrIntegral_sqrVs /Z0;
RFpulseProps.meanPower_W=RFpulseProps.energy_VAs /tau_pulse_s;
