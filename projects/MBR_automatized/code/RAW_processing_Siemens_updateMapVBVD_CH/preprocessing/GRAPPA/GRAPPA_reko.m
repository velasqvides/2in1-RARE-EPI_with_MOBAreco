function [Rk,coef] = GRAPPA_reko(Sk,ref,coef)
% function [Rk,coef] = GRAPPA_reko(Sk,ref,coef)
% Performs GRAPPA reconstruction of subsampled kSpace using a 3x2 kernel
%
% INPUT:
%       Sk: subsampled k-space
%       ref: reference scan
%       coefs: GRAPPA coefficients (if coef=0 then estimated)
%
% OUTPUT
%	Rk: Reconstructed k-space
%	coef: GRAPPA coefficients (weights)
%
% Implementation of
%
%     M. A. Griswold, P. M. Jakob, et. al..
%     Generalized autocalibrating partially parallel acquisitions
%     Mag. Reson. Med.   47(6):1202-1210.  Jun 2002
%
% Based on "recongrappa.m" by Scott Hoge  (shoge at ieee.org) and on code
% from Santiago Aja-Fernandez, LPI
%
% 2014 Matthias Dieringer
% matthias.dieringer@charite.de
%
% 2014-02-14, Till Huelnhagen, permuted the output image back to the input order
% 2014-02-14, Till Huelnhagen, avoided the use of function vec()

kernelX=[-1 0 1];

% --------------------------------
Sk=permute(Sk,[2,1,3]); %phase encoding data goes first
ref=permute(ref,[2,1,3]); %phase encoding data goes first
Mx=size(Sk,1);
My=size(Sk,2);
Ncoils=size(Sk,3);
Nl=size(ref,1); %number of ref lines
acs=(round(Mx/2-Nl/2+1)):round(Mx/2+Nl/2); %position of the reflines
Sk(acs,:,:)=ref(:,:,:); %put reflines into kspace
Slines=find(Sk(:,1,1)~=0);
% --------------------------------

%IF acs adjustent needed----
%if rem(rate,2)
%acs=acs(1:end-1);
%else
%acs=acs(2:end-1);
%end
%----------------------------

ImK = zeros( size(Sk,1), 1 );
ImK(Slines) = 1;  %Sampled lines set to 1

%ESTIMATE GRAPPA COEFFS --------------------------------------------
if coef==0
    QA = zeros(length(acs)*My, 6*Ncoils); %asumming 3x2 kernel
    coef = zeros(6*Ncoils,Ncoils);
    for ii = 1:length(acs),
        ind1 = acs(ii) + [-1 1];
        
        %% determine the k-space filling coefficients
        for jj=1:length(kernelX),
            y0 = mod( kernelX(jj) + [1:My],My);
            y0( y0 == 0 ) = My;
            QA( My*(ii-1) + (1:My), jj:length(kernelX):size(QA,2) ) = ...
                reshape( permute(squeeze( Sk(ind1,y0,:) ),[ 2 3 1 ]),My, Ncoils*length(ind1) );
        end;
    end;
    
    AA = QA'*QA;
    
    %Better with GNU cgsolv(), Copyright 2001 William Scott Hoge (shoge@ece.neu.edu or shoge@ieee.org)
    % taken from p.529 of Golub and Van Loan, "Matrix Computations," 3rd ed.
    % S. 10.2 The Conjugate Gradient Method
    % S. 10.2.6 Some Practical Details.
    %
    %n(:,l) = cgsolv( AA, A'*b, zeros(size(A,2),1), size(A,2) );
    %
    %SUBSTITUTION bicg
    
    
    for l=1:Ncoils
        b = squeeze(Sk( acs, :, l )).';
        b = b(:);
        [coef(:,l) ~]=bicg(AA,QA'*b);
        %TT Dummy variable to avoid warnings
    end
end %If COEF==0
%END ESTIMATION --------------------------------------------


%RECONSTRUCTION --------------------------------------------

for ii=1:length(ImK),  %For each line
    
    %1.- CONTROL OF INDEX
    if find(Slines==ii), continue, end
    
    %ELSE: if lines no sampled--> Recons
    tempo = [];
    ind1 = ii + [-1 1];
    if ( ind1(1) < 1) || (ind1(end) > Mx);
        klll=1;%do nothing
    else
        tempo = (ImK(ind1)==ones(length(ind1),1));
    end
    if ((min(ind1)<1)||(max(ind1)>Mx)||(sum(tempo)~=length(ind1)))
        continue;
    end;
    
    %2.- RECONSTRUCT K SPACE
    
    MA = zeros(My,3*2*Ncoils);
    for jj=1:3
        D1 = mod( kernelX(jj) + [1:My] , My );
        D1(D1==0) = My;
        MA((1:My),jj:3:3*2*Ncoils)=reshape(permute(squeeze(Sk(ind1,D1,:)),[ 2 3 1 ]),My, Ncoils*length(ind1));
    end %jj
    
    for l=1:Ncoils,
        Sk(ii,:,l)=MA*coef(:,l); %RECONST
    end
    
end  %ii
%Rk=Sk;
Rk=permute(Sk,[2,1,3]); %reestablish correct order of dimensions, TH, 2014-02-14

%END RECONSTR --------------------------------------------






