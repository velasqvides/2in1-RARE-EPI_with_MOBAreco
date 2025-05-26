function [kSpaceCorrected, kSpaceShifts] = correctGradientDelays(kSpace,gdelcalib,protPara)
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here
kSpaceCorrected = zeros(size(kSpace));
nSlices = protPara.nSlices;
kSpaceShifts = zeros(protPara.ETLfull,4,size(gdelcalib.gdelcalib_0,4),nSlices);

for slice = 1:nSlices
    if (protPara.ucTrajectory == 2) % 2 -> radial
        [kSpaceCorrected(:,:,:,slice), kSpaceShifts(:,:,:,slice)] = ...
            correctDelays(kSpace(:,:,:,slice),...
            gdelcalib,protPara,slice);
    end
end

    function [kSpaceCorrected, kSpaceShifts] = correctDelays(kSpaceUncorrected,gdelcalib,protPara,slice)
        % Function to correct gradient delays by estimating k-sapce shift with
        % spokes of opposite read direction accordint to the method described in
        % Block KT,Uecker M, Proc. 19th Ann. Meeting of ISMRM, Montr?al,CA,2011 (2816)

        reflect = true;

        nRARE = protPara.ETL_RARE;
        nEPI = protPara.ETL_EPI;
        ETL = protPara.ETLfull;

        if (protPara.interleaved == 1)
            interleaved = false;
            reversed = false;
        elseif (protPara.interleaved == 2)
            interleaved = true;
            reversed = false;
        elseif (protPara.interleaved == 3)
            interleaved = true;
            reversed = true;
        end

        % for CATFactor=0.1 this has to be changed so that nRARE=0
        % CATFactor  = 0.01;
        % ETL = MrProt.FastImaging.TurboFactor;
        % nRARE = floor(ETL*(CATFactor));
        if (size(gdelcalib.gdelcalib_0,5)>2)
            gdc0 =  permute( squeeze(gdelcalib.gdelcalib_0(:,:,:,:,slice)),[1 3 2]);
            gdc180 = permute( squeeze(gdelcalib.gdelcalib_180(:,:,:,:,slice)),[1 3 2]);
            gdc90 = permute( squeeze(gdelcalib.gdelcalib_90(:,:,:,:,slice)),[1 3 2]);
            gdc270 = permute( squeeze(gdelcalib.gdelcalib_270(:,:,:,:,slice)),[1 3 2]);
        else
            gdc0 = permute( gdelcalib.gdelcalib_0(:,:,:,:,end), [1 3 2]);
            gdc180 = permute( gdelcalib.gdelcalib_180(:,:,:,:,end), [1 3 2]);
            gdc90 = permute( gdelcalib.gdelcalib_90(:,:,:,:,end), [1 3 2]);
            gdc270 = permute( gdelcalib.gdelcalib_270(:,:,:,:,end), [1 3 2]);
        end

        [nFE,nSpokesFull,nCh]=size(kSpaceUncorrected);


        if (interleaved==true)
            % old version
            %index_EPI = (nRARE*(nSpokesFull/ETL)+1)+120:nSpokesFull;
            % corrected 06/30/2020
            index_EPI = (nRARE*(nSpokesFull/ETL)+1):nSpokesFull;
        else
            index_EPI = repmat(((nRARE+1):ETL)',1,nSpokesFull/ETL) + repmat(0:ETL:(nSpokesFull-ETL),ETL-nRARE,1);
            index_EPI = index_EPI(:);
        end

        if (interleaved && reversed)
            if ( mod(ETL-nRARE,2)==0 )
                reflectStart_img = 2;
            else
                reflectStart_img = 1;
            end
        else
            reflectStart_img = 1;
        end

        reflectStart_calib = 1;

        % for CATFactor=0.1 this has to be changed so that reflectStart is 2
        % reflectStart_calib = 2;
        % reflectStart_img = 2;

        nLinesGDelCalib = size(gdc0,2);
        index_EPI_GDelCalib = repmat(((nRARE+1):ETL)',1,nSpokesFull/ETL) + repmat(0:ETL:(nSpokesFull-ETL),ETL-nRARE,1);
        index_EPI_GDelCalib = index_EPI_GDelCalib(:);
        index_RARE_GDelCalib = repmat((1:nRARE)',1,nSpokesFull/ETL) + repmat(0:ETL:(nSpokesFull-ETL),nRARE,1);
        index_RARE_GDelCalib = index_RARE_GDelCalib(:);

        kSpaceUncorrected_reflected = kSpaceUncorrected;
        reflectIndexImg = index_EPI(reflectStart_img:2:end);
        reflectStart_img=1;
        if (mod(ETL-nRARE,2)~=0)
            tmp = zeros(nEPI,1);
            tmp(1:2:end) = 1;
            tmp = repmat(tmp, protPara.nSpokes,1);
            reflectIndexImg = index_EPI(tmp==1);
            %     tmpOdd = 1:2:size(index_EPI,2);
            %     tmpEven = 0:2:size(index_EPI,2);
            %     idx = repmat(repelem([false,true],ceil(nEPI/2)),1,size(tmpOdd,2)/(nEPI+1));
            %     tmpOdd(idx)= tmpEven(idx);
            %     %reflectIndexImg = index_EPI(tmpOdd);
            %     reflectIndexImg = index_EPI;
            %     clear tmpOdd tmpEven idx
        end
        if (reflect == true)
            kSpaceUncorrected_reflected(:,reflectIndexImg,:) = flipdim(kSpaceUncorrected(:,reflectIndexImg,:),1);
        end
        reflectIndexGDC = reflectStart_calib:2:floor((ETL-nRARE)*nLinesGDelCalib/ETL + .5);
        if (mod(ETL-nRARE,2)~=0)
            tmp = repelem(0:(nLinesGDelCalib/ETL-1),(ETL-nRARE+1)/2);
            reflectIndexGDC = (reflectStart_calib:2:2*size(tmp,2)) - tmp;
        end
        index_gdelcalib = index_EPI_GDelCalib(reflectIndexGDC);
        gdc0_reflected = gdc0;
        gdc90_reflected = gdc90;
        gdc180_reflected = gdc180;
        gdc270_reflected = gdc270;
        if (reflect == true)
            gdc0_reflected(:,index_gdelcalib,:)  = flipdim(gdc0_reflected(:,index_gdelcalib,:),1);
            gdc90_reflected(:,index_gdelcalib,:) = flipdim(gdc90_reflected(:,index_gdelcalib,:),1);
            gdc180_reflected(:,index_gdelcalib,:)= flipdim(gdc180_reflected(:,index_gdelcalib,:),1);
            gdc270_reflected(:,index_gdelcalib,:)= flipdim(gdc270_reflected(:,index_gdelcalib,:),1);
        end

        k_shifts_x=zeros(nCh,nLinesGDelCalib);
        k_shifts_y=zeros(nCh,nLinesGDelCalib);
        norm_spokes_x=zeros(nCh,nLinesGDelCalib);
        norm_spokes_y=zeros(nCh,nLinesGDelCalib);

        s_0 = abs(gdc0_reflected);
        s_180 = flipdim(abs(gdc180_reflected),1);
        s_90 = abs(gdc90_reflected);
        s_270 = flipdim(abs(gdc270_reflected),1);

        for Line=1:nLinesGDelCalib % performing correction for repeated acquisitions for averaging
            for Ch=1:nCh % performing correction for every channel
                % reflect 180? spoke and calculate magnitudes
                s_0_Ch_Line = s_0(:,Line,Ch);
                s_180_Ch_Line = s_180(:,Line,Ch);
                s_90_Ch_Line = s_90(:,Line,Ch);
                s_270_Ch_Line = s_270(:,Line,Ch);

                % fourier transform signal of spokes
                fft_s_0 = fftshift(ifft(ifftshift(s_0_Ch_Line),[],1));
                fft_s_180 = fftshift(ifft(ifftshift(s_180_Ch_Line),[],1));
                fft_s_90 = fftshift(ifft(ifftshift(s_90_Ch_Line),[],1));
                fft_s_270 = fftshift(ifft(ifftshift(s_270_Ch_Line),[],1));

                % calculate cross-correlation
                g_x=conj(fft_s_0).*(fft_s_180);
                g_y=conj(fft_s_90).*(fft_s_270);

                % linear regression of phase of g find support for linear
                % regression by signal drop of *threshold*% from maximum
                threshold = .2;
                % calculate range for x direction
                mag_fft_s_0 = abs(fft_s_0);
                [max_fft_s_0,max_fft_s_0_index] = max(mag_fft_s_0);
                range_lower_x = find(mag_fft_s_0(1:max_fft_s_0_index-1)<threshold*max_fft_s_0, 1, 'last' );
                range_upper_x = max_fft_s_0_index + find(fft_s_0(max_fft_s_0_index+1:end)<threshold*max_fft_s_0,1,'first');
                range_x = range_lower_x:range_upper_x;
                % calculate range for y direction
                mag_fft_s_90 = abs(fft_s_90);
                [max_fft_s_90,max_fft_s_90_index] = max(mag_fft_s_90);
                range_lower_y = find(mag_fft_s_90(1:max_fft_s_90_index-1)<threshold*max_fft_s_90, 1, 'last' );
                range_upper_y = max_fft_s_90_index + find(fft_s_90(max_fft_s_90_index+1:end)<threshold*max_fft_s_90,1,'first');
                range_y = range_lower_y:range_upper_y;

                % linear regression of g
                lin_regress_x = polyfit(range_x',angle(g_x(range_x,1)),1);
                lin_regress_y = polyfit(range_y',angle(g_y(range_y,1)),1);
                slope_x = lin_regress_x(1);
                slope_y = lin_regress_y(1);
                k_shifts_x(Ch,Line) = -slope_x*nFE/(4*pi);
                k_shifts_y(Ch,Line) = -slope_y*nFE/(4*pi);
                % plot phase in range used for fit (range depends on threshold)

                %figure();
                %plot(range_y,angle(g_y(range_y,1)),'blue'); hold on
                %plot(range_x,angle(g_x(range_x,1)),'red'); hold off

                % L2-Norm of 0? spokes for weighting to suppress influence of elements that recieve only
                % weak or no signal intensity
                norm_spokes_x(Ch,Line) = sqrt(sum(abs(s_0(:,Line,Ch)).^2));
                norm_spokes_y(Ch,Line) = sqrt(sum(abs(s_90(:,Line,Ch)).^2));
            end
        end

        % weighting with L2-Norm of spokes to suppress influence of elements that
        % recieve only weak or no signal intensity
        k_shifts_x_Ch_weighted = sum(k_shifts_x.*norm_spokes_x,1)./sum(norm_spokes_x,1);
        k_shifts_y_Ch_weighted = sum(k_shifts_y.*norm_spokes_y,1)./sum(norm_spokes_y,1);

        % average over Acquisitions
        %NAverages = wipPara.lNGdcScans; %nLinesGDelCalib/ETL
        NAverages = nLinesGDelCalib/ETL;
        k_shift_x = mean(k_shifts_x_Ch_weighted);
        k_shift_y = mean(k_shifts_y_Ch_weighted);
        k_shift_x_RARE = mean(k_shifts_x_Ch_weighted(index_RARE_GDelCalib(1:nRARE*NAverages)));
        k_shift_x_EPI = mean(k_shifts_x_Ch_weighted(index_EPI_GDelCalib(1:(ETL-nRARE)*NAverages)));
        k_shift_y_RARE = mean(k_shifts_y_Ch_weighted(index_RARE_GDelCalib(1:nRARE*NAverages)));
        k_shift_y_EPI = mean(k_shifts_y_Ch_weighted(index_EPI_GDelCalib(1:(ETL-nRARE)*NAverages)));

        %k_shift_x_EPI = k_shift_x_RARE;
        %k_shift_y_EPI = k_shift_y_RARE;

        k_shift_x_std = std(k_shifts_x_Ch_weighted);
        k_shift_y_std = std(k_shifts_y_Ch_weighted);
        k_shift_x_RARE_std = std(k_shifts_x_Ch_weighted(index_RARE_GDelCalib(1:nRARE*NAverages)));
        k_shift_x_EPI_std = std(k_shifts_x_Ch_weighted(index_EPI_GDelCalib(1:(ETL-nRARE)*NAverages)));
        k_shift_y_RARE_std = std(k_shifts_y_Ch_weighted(index_RARE_GDelCalib(1:nRARE*NAverages)));
        k_shift_y_EPI_std = std(k_shifts_y_Ch_weighted(index_EPI_GDelCalib(1:(ETL-nRARE)*NAverages)));


        k_shift_x_ETLLine_mean = mean(reshape(k_shifts_x_Ch_weighted,ETL,NAverages),2);
        k_shift_y_ETLLine_mean = mean(reshape(k_shifts_y_Ch_weighted,ETL,NAverages),2);
        k_shift_x_ETLLine_sd = std(reshape(k_shifts_x_Ch_weighted,ETL,NAverages),[],2);
        k_shift_y_ETLLine_sd = std(reshape(k_shifts_y_Ch_weighted,ETL,NAverages),[],2);

        if (interleaved && reversed)
            echoIndex=zeros(1,ETL);
            for i=1:ETL/2
                echoIndex(1,2*i-1) = i;
                echoIndex(1,2*i) = ETL-i+1;
            end
        else
            echoIndex=1:ETL;
        end

        %k_shifts_x_mean = k_shift_x_ETLLine_mean(echoIndex);
        %k_shifts_y_mean = k_shift_y_ETLLine_mean(echoIndex);
        % use k_shift averaged over echo train
        %k_shifts_x_mean = repelem(mean(k_shift_x_ETLLine_mean(echoIndex)),ETL)';
        %k_shifts_y_mean = repelem(mean(k_shift_y_ETLLine_mean(echoIndex)),ETL)';
        % use averaged k_shift of RARE echoes
        k_shifts_x_mean = repelem(mean(k_shift_x_ETLLine_mean(1:nRARE)),ETL)';
        k_shifts_y_mean = repelem(mean(k_shift_y_ETLLine_mean(1:nRARE)),ETL)';

        % use averaged k_shifts of RARE/EPI echoes
        k_shifts_x_mean_RARE = repelem(mean(k_shift_x_ETLLine_mean(1:nRARE)),ETL)';
        k_shifts_y_mean_RARE = repelem(mean(k_shift_y_ETLLine_mean(1:nRARE)),ETL)';
        k_shifts_x_mean_EPI  = repelem(mean(k_shift_x_ETLLine_mean(1:nRARE)),ETL)';
        k_shifts_y_mean_EPI  = repelem(mean(k_shift_y_ETLLine_mean(1:nRARE)),ETL)';

        %--------------------------------------------------------------------------
        % correction with phasefactor pr echo in echo train
        %--------------------------------------------------------------------------
        %directly use avergaed k-space shifts
        k_shifts_x_mean_array_RARE = permute( repmat(k_shifts_x_mean_RARE(1:nRARE),[nSpokesFull/ETL 1 nFE]), [3 1 2] );
        k_shifts_y_mean_array_RARE = permute( repmat(k_shifts_y_mean_RARE(1:nRARE),[nSpokesFull/ETL 1 nFE]), [3 1 2] );
        k_shifts_x_mean_array_EPI = permute( repmat(k_shifts_x_mean_EPI(1:nEPI),[nSpokesFull/ETL 1 nFE]), [3 1 2] );
        k_shifts_y_mean_array_EPI = permute( repmat(k_shifts_y_mean_EPI(1:nEPI),[nSpokesFull/ETL 1 nFE]), [3 1 2] );
        %directly use k-space shifts from each echo
        % k_shifts_x_mean_array_RARE = permute( repmat(k_shift_x_ETLLine_mean(1:nRARE),[nSpokesFull/ETL 1 nFE]), [3 1 2] );
        % k_shifts_y_mean_array_RARE = permute( repmat(k_shift_y_ETLLine_mean(1:nRARE),[nSpokesFull/ETL 1 nFE]), [3 1 2] );
        % k_shifts_x_mean_array_EPI = permute( repmat(k_shift_x_ETLLine_mean((nRARE+1):ETL),[nSpokesFull/ETL 1 nFE]), [3 1 2] );
        % k_shifts_y_mean_array_EPI = permute( repmat(k_shift_y_ETLLine_mean((nRARE+1):ETL),[nSpokesFull/ETL 1 nFE]), [3 1 2] );
        if (protPara.bSameSpokeInET==1) % 1 -> false, 2 -> true
            angles_RARE = repmat((0:nSpokesFull/ETL*nRARE-1)*pi/(nSpokesFull/ETL*nRARE),[nFE 1 nCh]);
            angles_EPI = repmat((0:nSpokesFull/ETL*nEPI-1)*pi/(nSpokesFull/ETL*nEPI),[nFE 1 nCh]);
        else % case of no blip gradients
            angles = repmat(reshape(repmat(0:(nSpokesFull/ETL-1),ETL,1),1,nSpokesFull)*pi/(nSpokesFull/ETL),[nFE 1 nCh]);
                end
                % interpolate k-space shift for angle between 0degree 90degree
                if (protPara.bSameSpokeInET==1) % 1 -> false, 2 -> true
                    k_shifts_RARE = ((cos(2*angles_RARE)+1).*(k_shifts_x_mean_array_RARE+k_shifts_x_mean_array_RARE)/2+ + (-cos(2*angles_RARE)+1).*(k_shifts_y_mean_array_RARE+k_shifts_y_mean_array_RARE)/2 )/2;
                    k_shifts_EPI = ((cos(2*angles_EPI)+1).*(k_shifts_x_mean_array_EPI+k_shifts_x_mean_array_EPI)/2 + (-cos(2*angles_EPI)+1).*(k_shifts_y_mean_array_EPI+k_shifts_y_mean_array_EPI)/2 )/2;
                    k_shifts = [k_shifts_RARE k_shifts_EPI]; %#ok<NASGU>
                else
                    k_shifts = ((cos(2*angles)+1).*k_shifts_x_mean_array_EPI + (-cos(2*angles)+1).*k_shifts_y_mean_array_EPI )/2;
                    %k_shifts = ((cos(2*angles)+1).*(k_shifts_x_mean_array+k_shifts_y_mean_array)/2 + (-cos(2*angles)+1).*(k_shifts_x_mean_array+k_shifts_y_mean_array)/2 )/2;
                end
                phasefactor = exp(-sqrt(-1)*2*pi/nFE*k_shifts.*repmat((1:nFE)',[1 nSpokesFull nCh]));
                %--------------------------------------------------------------------------

                % old correction not per echo
                %phi = (0:(nSpokesFull-1))*pi/nSpokesFull;
                %k_shift_phi_RARE_EPI = ((cos(2*phi)+1)*k_shift_x_RARE + (-cos(2*phi)+1)*k_shift_y_RARE )/2;
                %k_shift_phi_RARE_EPI(1,index_EPI) = ((cos(2*phi(1,index_EPI))+1)*k_shift_x_EPI + (-cos(2*phi(1,index_EPI))+1)*k_shift_y_EPI )/2;
                %k_shift_phi_RARE_EPI_array = reshape(repmat(k_shift_phi_RARE_EPI,nFE,nCh),nFE,nSpokesFull,nCh);
                %phasefactor_matrix = reshape(repmat(1:nFE,nSpokesFull*nCh,1)',nFE,nSpokesFull,nCh);
                %phasefactor = exp(-sqrt(-1)*2*pi/nFE*k_shift_phi_RARE_EPI_array.*phasefactor_matrix);

                % shift spokes in k-space
                fft_spokes_reflected = fftshift(ifft(ifftshift((kSpaceUncorrected_reflected)),[],1));
                kSpaceCorrected = ifftshift(fft(fftshift(fft_spokes_reflected.*phasefactor),[],1));
                if (reflect==true)
                    kSpaceCorrected(:,reflectIndexImg,:) = flipdim(kSpaceCorrected(:,reflectIndexImg,:),1);
                    %kSpaceCorrected(:,index_EPI(reflectStart_img:2:floor(nSpokesFull*(1-CATFactor)+0.5)),:) = flipdim(kSpaceCorrected(:,index_EPI(reflectStart_img:2:floor(nSpokesFull*(1-CATFactor)+0.5)),:),1);
                end


                kSpaceShifts = [k_shift_x_ETLLine_mean k_shift_x_ETLLine_sd k_shift_y_ETLLine_mean k_shift_y_ETLLine_sd];

    end % end nested function
end % end main function