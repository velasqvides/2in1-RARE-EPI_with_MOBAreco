%call function with: [temp_pulse_playedout_V,temp_PulsProps]=calculatePulseProperties_fromSimulatedPulse(RealAmp_1,filename_1stRFpulse);
function RFpulseProps=calculatePulseProperties_fromSimulatedPulse(filename_1stRFpulse,RealAmp_1)

%%
A=importdata(filename_1stRFpulse);
pulse_playedout_V=A.data(logical(A.data)); % select only data unequal 0
dt=1e-6; % 1us
tau_pulse_s=length(pulse_playedout_V)*dt;
amp_pulse_V=max(pulse_playedout_V);

if RealAmp_1~=amp_pulse_V  
    disp('amplitudes from pulse file does not match the amplitude from the protocol:')
    disp(['protocol amplitude: ',num2str(RealAmp_1),' V'])
    disp(['pulse file amplitude: ',num2str(amp_pulse_V),' V'])
    disp('used amplitude for calculation: protocol amplitude')
    pulse_playedout_V=pulse_playedout_V/amp_pulse_V*RealAmp_1;
    amp_pulse_V=RealAmp_1;
end

%%
Z0=50;%Ohm=V/A; impedance

RFpulseProps.sourceOfRFpusle='RF pulse simulated with POET';
RFpulseProps.filenameOfRFpulse=filename_1stRFpulse;

RFpulseProps.peakAmplitude_V=amp_pulse_V;
RFpulseProps.pulse_duration_s=tau_pulse_s;

RFpulseProps.pulse_playedout_V = pulse_playedout_V;%[V]
RFpulseProps.rastertime_s=dt;

RFpulseProps.pulseIntegral_Vs = sum(pulse_playedout_V.*dt);%unit [V*s]
RFpulseProps.pulseSqrIntegral_sqrVs = sum(pulse_playedout_V.^2 .*dt);%unit [V^2*s]
RFpulseProps.energy_VAs=RFpulseProps.pulseSqrIntegral_sqrVs /Z0;
RFpulseProps.meanPower_W=RFpulseProps.energy_VAs /tau_pulse_s;

