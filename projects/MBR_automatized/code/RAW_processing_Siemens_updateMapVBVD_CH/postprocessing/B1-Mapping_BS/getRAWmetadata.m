function [Pulse,Sequence]=getRAWmetadata(MrProt)
% This function gets information about the protocol and the RF pulses used
% for B1-Mapping. Some of the parameters are not used in the calculation but
% are rather for just your information.
%
% INPUT: Siemens RAW data file including path.
% OUTPUT: PULSE and SEQUENCE details from the raw data header
% 
% 06-2013 Matthias Dieringer
% last updated 01-2014
% matthias.dieringer@charite.de


% PULSE PARAMETERS
if MrProt.WiPMemBlock.alFree(1)==0
    Pulse.TYPE='rect';
elseif MrProt.WiPMemBlock.alFree(1)==1
    Pulse.TYPE='fermi';
elseif MrProt.WiPMemBlock.alFree(1)==2
    Pulse.TYPE='adiabatic';
end

Pulse.TotalDuration=MrProt.WiPMemBlock.alFree(3); %duration of the B1 pulse
Pulse.SEG=MrProt.WiPMemBlock.alFree(2); %segments (samples)of the B1 pulse
Pulse.OLD=MrProt.WiPMemBlock.adFree(1); % single loop duration
Pulse.HL=MrProt.WiPMemBlock.adFree(2); %half loops of the B1 pulse
Pulse.P0=MrProt.WiPMemBlock.adFree(3); % initial excitation pulse flip angle
Pulse.PF=MrProt.WiPMemBlock.adFree(4); % flip back pulse flip angle
Pulse.LoopFA=MrProt.WiPMemBlock.adFree(5); % B1 pulse flip angle
if(exist ('MrProt.WiPMemBlock.adFree(6)'))
    Pulse.PS2D=MrProt.WiPMemBlock.adFree(6); % sinc RF pulse phase advance
else
    Pulse.PS2D=0;
end
if(exist('MrProt.WiPMemBlock.adFree(7)'))
    Pulse.PS3D=MrProt.WiPMemBlock.adFree(7); % excitation pulse phase advance
else
    Pulse.PS3D=0;
end

Pulse.TotalFA=Pulse.P0+Pulse.LoopFA+Pulse.PF; %total flip angle
Pulse.NominalB1=1E6*(Pulse.LoopFA*pi/180)/(Pulse.TotalDuration*1E-6*2.67522212E8);

% OTHER PARAMETERS USED IN THE PROTOCOL
Sequence.Nucleus=MrProt.TXSPEC.NucleusInfo(1).Nucleus;
Sequence.Slices=MrProt.SliceArray.Size;
Sequence.FOV(1:Sequence.Slices,1)=MrProt.SliceArray.Slice.PhaseFOV;
Sequence.Thickness(1:Sequence.Slices,1)=MrProt.SliceArray.Slice.Thickness;
Sequence.Partitions=MrProt.KSpace.Partitions;
Sequence.Resolution=MrProt.KSpace.BaseResolution;
Sequence.PhaseEncLines=MrProt.KSpace.PhaseEncodingLines;
Sequence.Averages=MrProt.Averages;
Sequence.Contrasts=MrProt.Contrasts;
Sequence.TR=MrProt.TR*1E-3;
Sequence.TA=MrProt.TotalScanTimeSec;
Sequence.TE1=MrProt.TE(1)*1E-6;
Sequence.TE2=MrProt.TE(2)*1E-6;
Sequence.DeltaTE=Sequence.TE2-Sequence.TE1;

% REFERENCE VOLTAGE OF THE FLIP ANGLE
Sequence.RefVolt(1,:)=MrProt.TXSPEC.NucleusInfo(1).ReferenceAmplitude;
Sequence.RealVolt(1,:)=MrProt.TXSPEC.RFPULSE.Amplitude;
Sequence.ReqVolt(1,:)=Sequence.RefVolt*(2*Pulse.LoopFA/Pulse.HL/180)*(1E+3/Pulse.OLD);
Sequence.FA=MrProt.FlipAngleDegree;
Sequence.ApproxFA=Sequence.FA*Sequence.RealVolt/Sequence.ReqVolt;

% use initial excitation pulse as p0 if this is not the COMPO method
% (for clarity, we should rather use an extra parameter)
if ~(strcmp(Pulse.TYPE,'rect'))
    Pulse.P0=Sequence.FA;
end

% for 2D we have an extra sinc pulse and no flip back pulse
if (MrProt.KSpace.Dimension==2)
    Sequence.Partitions=Sequence.Slices;
    Pulse.PF=0;
    Pulse.TotalFA=Pulse.P0+Pulse.LoopFA;
end

% 2D image imterpolation (zero filling) yes/no?
Sequence.ImageInterpolation=isfield(MrProt.KSpace,'Interpolation'); 
end