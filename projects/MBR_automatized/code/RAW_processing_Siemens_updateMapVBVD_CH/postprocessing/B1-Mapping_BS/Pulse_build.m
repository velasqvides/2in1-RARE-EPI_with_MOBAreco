function [dur,CmplxPulse,CmplxPulseConj]=Pulse_build(Pulse)
% function Pulse_build(Pulse)
%
% INPUT:
% Pulse=RF pulse struct containing B1 pulse details
%
% OUTPUT:
% dur=pulse duration [us]
% CmplxPulse=Complex RF pulse (magn and phase) [Gauss]
% CmplxPulseConj=Conjugate complex RF pulse (magn and phase) [Gauss]
%
% simplified, vectorized, and adapted to incorporate adiabatic Bloch
% Siegert pulse, 2013 Matthias Dieringer
% last updated 01-2014
% modified 2016-09-21: now capable for X-Nuclei (Olli Weinberger)

% generate complex RF pulse envelops dependend on the selected method
if strcmp(Pulse.TYPE,'adiabatic') % adiabatic Bloch-Siegert
    load('ABS.mat')
    CmplxPulse=conj(ABS);
    % calculate duration of each time step dependend on the pulse duration
    dur=ones(numel(CmplxPulse),1)*Pulse.TotalDuration/numel(CmplxPulse);
else
    if strcmp(Pulse.TYPE,'rect') %hard pulse (PhiFA-Cup composite RF pulse)
        flip=ones(1,Pulse.TotalDuration);
        % n is a vector centered at zero ranging from -0.5 to almost 0.5 in
        % steps of Pulse.SEG and with number of samples=Pulse.TotalDuration
        % (for the fermi pulse this increase is continous)
        n=floor(Pulse.SEG*((1:Pulse.TotalDuration)/Pulse.TotalDuration))/Pulse.SEG-0.5;
    elseif strcmp(Pulse.TYPE,'fermi')
        % pulse parameters FERMI
        fSteep=10.0;
        fSigma=0.3;
        % monotonic increasing vector n centered at zero
        n=((0:Pulse.TotalDuration-1)/Pulse.TotalDuration)-0.5;
        flip=1./(1.0+exp(fSteep*(abs(n)-fSigma)/fSigma));
    end
    
    % calculate phase course of the RF pulse (this formula is according to
    % the formula used in the sequence)
    pha = 2*pi*((n+.5)*(Pulse.HL/2)+.5+Pulse.PS3D/360-.5);
% Phi-FA:
%     pha= (360*((n+.5)*(Pulse.HL/2)+.5+Pulse.PS3D/360-.5)-157.5)/180*pi;
    
    CmplxPulse=(flip.*exp(1i*pha))';%now: max(CmplxPulse)=1
    % calculate duration of each time step
    dur=ones(numel(CmplxPulse),1);
end

%calculate integral of the B1 pulse
loopINT=sum(abs(CmplxPulse));%flip(1+samplesP0:end));

% normalize pulse amplitude to integral and scale it to flip angle
CmplxPulse=CmplxPulse/loopINT*Pulse.LoopFA;%now: CmplxPulse has the unit degree

% This is the initial excitation pulse. In the PhiFA-Cup method this is
% incorporated in the hard pulse. In 2D Methods this is the sinc RF pulse.
% We put this pulse in front of the B1 pulse with -90? phase shift.
if Pulse.P0~=0
    %simply assume 500us RF pulse, this is arbitrary but non-critical
    dur=[500; dur];
    CmplxPulse=[Pulse.P0*exp(-1i*pi/2); CmplxPulse];
%     CmplxPulse=[Pulse.P0; CmplxPulse];
end

% If a flip back pulse exists, then add it to the back... (has +90? phase
% shift!)
if Pulse.PF
    %assume 500us RF pulse
    dur=[dur; 500];
    CmplxPulse=[CmplxPulse; Pulse.PF*exp(1i*pi/2)];
end

% until here the integrals represent the flip angle in degrees, now convert to
% Gauss for bloch simulation

gamma=Pulse.gamma_Nucleus_rad_s_T;

%CmplxPulse=1E12*CmplxPulse*pi./(180*gamma*dur)*1E-2; %[uGauss]

CmplxPulse=CmplxPulse*pi/180./dur/gamma * 1E10; %[Gauss]

%complex conjugate of pulse in Gauss (only conjugate B1 relevant pulse)
CmplxPulseConj=CmplxPulse;
% CmplxPulseConj(1+(Pulse.P0~=0):end-(Pulse.PF~=0))=conj(CmplxPulse(1+(Puls
% e.P0~=0):end-(Pulse.PF~=0)));
CmplxPulseConj=conj(CmplxPulse);

end
