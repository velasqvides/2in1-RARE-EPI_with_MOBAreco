function  AFI = AFI_read_TXfct(varargin)
 %% ************************************************************************
%
%   Reads AFI maps from DCM images and returns processed B1+ maps
%   
%   Author: S.Schmitter
%   Date: April 2018
%
%
%   Dependencies:
%   - parseVariableInputs
%   - parse_Siemens_Shadow
%
%   Modified: line 83: alpha = rad2deg(abs(acos((r*n-1)./(n-r)))).*~excludemask;
%
%   INPUT:                                                  [unit]
%   ----------------------------------------------------------------
%   
%
%   Options (with standard prefs)                           [unit]
%   ----------------------------------------------------------------
%   optsdef.THRESHOLD     = 5;  cutoff threshold; 
%                               in multiples of noise level
%   optsdef.ALPHATHRESH   = 4;  lowest alpha value displayed in deg
%
%   OUTPUT:
%   ----------------------------------------------------------------
%   AFI     structure with map and options
%
%% ************************************************************************

    %default options
    optsdef.THRESHOLD       = 5; % in units of mean noise level
    optsdef.ALPHATHRESH     = 5; % in degrees

    opts    = catstruct(optsdef,parseVariableInputs(varargin));


    [sDcms,sDcmInfo] = readDcmFolder();
    %from eddie:
    [img, ser, mrprot] = parse_siemens_shadow_V2(sDcmInfo(1));

    %alpha = arccos((r*n-1)/(n-1))
    %r = S2/S1
    %n = TR2/RE1

    %get mean repetition time from dicominfo field:
    TRmean = double(sDcmInfo(1).RepetitionTime)*1E3;
    TRdiff = double(mrprot.sWiPMemBlock.alFree(12)); %mrprot field with difference

    TR1 = TRmean - TRdiff;
    TR2 = TRmean + TRdiff;

    n = TR2/TR1;

    r = double(sDcms(:,:,2:2:end))./double(sDcms(:,:,1:2:end));

    %calculate the mask:
    %-------------------
    %take images from TR2:
    dcm2 = double(double(sDcms(:,:,2:2:end)));

    dcm2cent = double(dcm2(:,:,end/2));
    %corners of the image
    xx = [2,3];
    yy = [2,3];

    ntmp(1) = mean(mean(dcm2cent(xx,yy)));
    ntmp(2) = mean(mean(dcm2cent(end-xx,yy)));
    ntmp(3) = mean(mean(dcm2cent(xx,end-yy)));
    ntmp(4) = mean(mean(dcm2cent(end-xx,end-yy)));

    %noise level:
    noise = mean(ntmp);

    %here is the mask
    excludemask = dcm2 < noise*opts.THRESHOLD;

    %claculate the map:
    %------------------

    %this is the b1map:
    alpha = rad2deg(abs(acos((r*n-1)./(n-r)))).*~excludemask;
    %mod 20180404
    disp('new AFI recon');
    alphamask = alpha < 4;
    alpha = alpha.*~alphamask;

    if(isfield(mrprot.sSliceArray.asSlice,'sPosition'))
        mpos = mrprot.sSliceArray.asSlice.sPosition;

        if(isfield(mpos,'dSag'))
            dSag = mpos.dSag;
        else
            dSag = 0;
        end
        if(isfield(mpos,'dTra'))
            dTra = mpos.dTra;
        else
            dTra = 0;
        end
        if(isfield(mpos,'dCor'))
            dCor = mpos.dCor;
        else
            dCor = 0;
        end
    else
        %modss 20160524
        dSag = 0;
        dTra = 0;
        dCor = 0;
    end
    
    %info structure
    Info{1}.name = 'Matrix size'; 
    Info{1}.value = num2str(size(alpha));
    Info{2}.name = 'TR1'; 
    Info{2}.value = [num2str(round(TR1/1E3)),' ms'];
    Info{3}.name = 'TR2'; 
    Info{3}.value = [num2str(round(TR2/1E3)),' ms'];
    Info{4}.name = 'n (ratio)'; 
    Info{4}.value = num2str(n);
    Info{5}.name = 'FA (nominell)'; 
    Info{5}.value = [num2str(mrprot.adFlipAngleDegree),' deg'];
    Info{6}.name = 'U_ref'; 
    Info{6}.value = [num2str(mrprot.sTXSPEC.asNucleusInfo(1).flReferenceAmplitude),' V'];
    Info{7}.name = 'Position (Tra/Sag/Cor)'; 
    Info{7}.value = [num2str(dTra),' ',num2str(dSag),' ',num2str(dCor),' mm'];
    Info{8}.name = 'Position (Siemens)'; 
    Info{8}.value = sDcmInfo(1).Private_0051_100e;
    
    
    AFI.InfoID = 'myinfo_AFI';
    
    AFI.alpha = alpha;
    AFI.noise = noise;
    AFI.opts = opts;
    AFI.TR1 = TR1;
    AFI.TR2 = TR2;
    AFI.n = n;
    AFI.sDcmInfo = sDcmInfo(1);
    AFI.mrprot = mrprot;
    AFI.Info = Info;
    AFI.DCM_TR1 = sDcms(:,:,1:2:end);
    AFI.DCM_TR2 = sDcms(:,:,2:2:end);

end
