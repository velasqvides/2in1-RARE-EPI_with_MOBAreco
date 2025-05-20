%call function with: [temp_pulse_playedout_V,temp_PulsProps]=calculatePulseProperties_fromMeasuredPulse(RealAmp_1,filename_1stRFpulse);
function RFpulseProps=calculatePulseProperties_fromMeasuredPulse(filename_1stRFpulse,RealAmp_1)

%%
load(filename_1stRFpulse);
pulseEnvelope=env;
dt=Tinterval;%in s

maxPulse=max(pulseEnvelope);
pulse_playedout_V=RealAmp_1/maxPulse .*pulseEnvelope;
tau_pulse_s=length(pulse_playedout_V)*dt;
%%
Z0=50;%Ohm=V/A; impedance

RFpulseProps.sourceOfRFpusle='measured RF pulse';
RFpulseProps.filenameOfRFpulse=filename_1stRFpulse;

RFpulseProps.peakAmplitude_V=RealAmp_1;
RFpulseProps.pulse_duration_s=tau_pulse_s;

RFpulseProps.pulse_playedout_V = pulse_playedout_V;%[V]
RFpulseProps.rastertime_s=dt;

RFpulseProps.pulseIntegral_Vs = sum(pulse_playedout_V.*dt);%unit [V*s]
RFpulseProps.pulseSqrIntegral_sqrVs = sum(pulse_playedout_V.^2 .*dt);%unit [V^2*s]
RFpulseProps.energy_VAs=RFpulseProps.pulseSqrIntegral_sqrVs /Z0;
RFpulseProps.meanPower_W=RFpulseProps.energy_VAs /tau_pulse_s;


