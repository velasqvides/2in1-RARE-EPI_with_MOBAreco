function Data=Sensitivity_build_mex(Pulse,Sequence,Data)
% function [Curve]=Sensitivity_build_mex(Pulse,Data)
%
% INPUT: Pulse = B1+ sensitive pulse details
%        Data = struct containing simulation details
% OUTPUT:Curve = struct containing interpolation curve and details
%
% Builds an interpolation curve based on the B1+ pulse setting in the
% sequence and simulates magnetization with positive and negative
% freq offset.
% 2013 Matthias Dieringer
% last updated 01-2014
% modified 2016-09-21: Now capable for X-Nuclei (Olli Weinberger)
% modified 2018-10-05: Add sequence info with nucleus info (Thomas Eigentler)


% build RF pulse traces

[dur_us,CmplxPulse,CmplxPulseConj]=Pulse_build(Pulse);

% define span of B0 and B1 to simulate
Data.B0vec=2*((0:Data.SpinsB0-1)/(Data.SpinsB0-1)-0.5)*Data.DistrB0;
%Data.B1vec=Data.B1Max*(0:Data.SpinsB1-1)/(Data.SpinsB1);
Data.B1vec=Data.B1Max*(1:Data.SpinsB1)/(Data.SpinsB1);
[Data.B1mat,Data.B0mat]=meshgrid(Data.B1vec,Data.B0vec);
dur=dur_us/1E6;%[s]
if isfield(Data, 'T1')
    T1=Data.T1;
else
    T1=200;%[s]
end

if isfield(Data, 'T2')
    T2=Data.T2;
else
    T2=200;%[s]
end

% Bloch simulation of magnetization using defined RF pulse traces
if strcmp(Sequence.Nucleus,'1H')                                            %TE: Sequence Info added                                   
    for u=1:Data.SpinsB1
        %positive freq offset
        %[mx,my,mz]                  = bloch(b1,                            gr,     tp,t1,t2,df,        dp,mode,mx,my,mz)
        [MX1(:,u),MY1(:,u),MZ1(:,u)] = bloch(Data.B1vec(u).*CmplxPulse    ,dur.*0,dur,T1,T2,Data.B0vec,1,0);
        %negative freq offset
        [MX2(:,u),MY2(:,u),MZ2(:,u)] = bloch(Data.B1vec(u).*CmplxPulseConj,dur.*0,dur,T1,T2,Data.B0vec,1,0);
    end
elseif strcmp(Sequence.Nucleus,'23Na')
    for u=1:Data.SpinsB1
        %positive freq offset
        %[mx,my,mz]                  = bloch_Na(b1,                            gr,     tp,t1,t2,df,        dp,mode,mx,my,mz)
        [MX1(:,u),MY1(:,u),MZ1(:,u)] = bloch_Na(Data.B1vec(u).*CmplxPulse    ,dur.*0,dur,T1,T2,Data.B0vec,1,0);
        %negative freq offset
        [MX2(:,u),MY2(:,u),MZ2(:,u)] = bloch_Na(Data.B1vec(u).*CmplxPulseConj,dur.*0,dur,T1,T2,Data.B0vec,1,0);
    end
end

% for u=1:Data.SpinsB1
%     %positive freq offset
%     %[mx,my,mz]                  = bloch(b1,                            gr,     tp,t1,t2,df,        dp,mode,mx,my,mz)
%     ##[MX1(:,u),MY1(:,u),MZ1(:,u)] = bloch(Curve.B1vec(u).*CmplxPulse    ,dur.*0,dur,T1,T2,Curve.B0vec,1,0);
%     %negative freq offset
%     [MX2(:,u),MY2(:,u),MZ2(:,u)] = bloch(Curve.B1vec(u).*CmplxPulseConj,dur.*0,dur,T1,T2,Curve.B0vec,1,0);
% end


Phase1=atan2(MX1,MY1);
Data.Magnitude1=sqrt(MY1.^2+MX1.^2);
Phase2=atan2(MX2,MY2);
Data.Magnitude2=sqrt(MY2.^2+MX2.^2);

Data.DeltaPhi=(Phase1-Phase2)*180/pi; %phase unwrap works with degrees here
Data.DeltaPhi=UnWrapMatBest(Data.DeltaPhi,size(Data.DeltaPhi,1)/2,1,200,100,1)+180;
end
