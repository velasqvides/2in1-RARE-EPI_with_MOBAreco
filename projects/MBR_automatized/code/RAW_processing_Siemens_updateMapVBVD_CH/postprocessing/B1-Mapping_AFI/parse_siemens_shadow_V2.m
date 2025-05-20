function [img, ser, mrprot] = parse_siemens_shadow_V2(dcm)
% [img, ser, mrprot] = parse_siemens_shadow(dcm)
% function to parse siemens numaris 4 shadow data
% returns three structs with image, series header, mrprot info
% does not work with arrayed dcm()
%    dependencies: parse_mrprot.m
%                  c_str.m
%                  mread.m

if (size(dcm,2) > 1)
    error('parse_siemens_shadow does not work on arrayed dicominfo data!')
end

ver_string = dcm.Private_0029_1008;
csa_string = dcm.Private_0029_10xx_Creator;

%mod ssm 20161006
if (strcmp(ver_string,'IMAGE NUM 4') || strcmp(char(ver_string(1:end-1).'),'IMAGE NUM 4') || strcmp(char(ver_string.'),'IMAGE NUM 4'))
%if (strcmp(char(ver_string(1:end-1).'),'IMAGE NUM 4'))
    if (strcmp(csa_string,'SIEMENS CSA HEADER'))
        img = parse_shadow_func(dcm.Private_0029_1010);
        ser = parse_shadow_func(dcm.Private_0029_1020);
    else
        error('shadow: Invalid CSA HEADER identifier: %s',csa_string);
    end
elseif (strcmp(ver_string,'SPEC NUM 4'))
    if (strcmp(csa_string,'SIEMENS CSA NON-IMAGE'))
        if isfield(dcm,'Private_0029_1210')
            img = parse_shadow_func(dcm.Private_0029_1210);
            ser = parse_shadow_func(dcm.Private_0029_1220);
        else %VB13
            img = parse_shadow_func(dcm.Private_0029_1110);
            ser = parse_shadow_func(dcm.Private_0029_1120);
        end
    else
        error('shadow: Invalid CSA HEADER identifier: %s',csa_string);
    end
else
    error('shadow: Unknown/invalid NUMARIS version: %s',ver_string);
end

% now parse the mrprotocol
if isfield(ser, 'MrPhoenixProtocol') % VB13
    MrProtocol = char(ser.MrPhoenixProtocol);
else
    MrProtocol = char(ser.MrProtocol);
end
spos = strfind(MrProtocol,'### ASCCONV BEGIN ###');
epos = strfind(MrProtocol,'### ASCCONV END ###');
%mod ssm for VD13D 
%original: only else expression
if(isempty(spos))
    spos = strfind(MrProtocol,'### ASCCONV BEGIN');
    tmp = MrProtocol(spos+18:end);
    spos = strfind(tmp,'###');
    epos = strfind(tmp,'### ASCCONV END ###');
    MrProtocol = tmp(spos+4:epos-2);
else
    MrProtocol = MrProtocol(spos+22:epos-2);
end
mrprot = parse_mrprot(MrProtocol);

%--------------------------------------------------------------------------

function hdr = parse_shadow_func(dcm)
% internal function to parse shadow header

% input (dcm) is uint8 using little endian ordering; since this could be
% run on a little or big endian machine, we need to interpret

% scan through the data byte by byte

%fprintf('\nReading IMAGE header\n');
fp = 9;                                                                 % skip 4 chars 'SV10' + unknown int32
[nelem, fp] = mread(dcm, fp, 1, 'int32-le');                            % # of elements int32
fp = fp + 4;                                                            % skip unknown int32 (77)
%fprintf('Found %d elements\n', nelem);
for y=1:nelem
    %data_start_pos = fp;
    [tag, fp] = mread(dcm, fp, 64, 'c_str');                            % field name c_str[64]
    tag = strrep(tag,'-','_'); % remove invalid chars from field name
    [vm, fp] = mread(dcm, fp, 1, 'int32-le');
    [vr, fp] = mread(dcm, fp, 4, 'c_str');                              % vr string c_str[4]
    fp = fp + 4;                                                        % skip SyngoDT int32
    [NoOfItems, fp] = mread(dcm, fp, 1, 'int32-le');                    % NoOfItems int32

    if (vm == 0) % this can happen in spectroscopy files, VB13 image files
        vm = NoOfItems;
    end

    %str_data = '';

    if (NoOfItems > 1)
        fp = fp + 4;                                                    % skip unknown int32 (77)

        for z=1:vm
            [int_items, fp] = mread(dcm, fp, 4, 'int32-le');            % field width
            use_fieldwidth = ceil(int_items(4) / 4) * 4;
            [tmp_data, fp] = mread(dcm, fp, use_fieldwidth, 'c_str');
            %str_data = [str_data tmp_data];
            %if (z < vm), str_data = [str_data '\']; end

            switch vr
                case {'AE','AS','CS','DA','DT','LO','LT','OB','OW','PN','SH','SQ','ST','TM','UI','UN','UT'}
                    % these are string values
                    %fprintf('String VR %s, data = %s\n',vr,str_data);
                    if (z == 1), val_data = cell(vm,1); end
                    val_data{z} = tmp_data;
                case {'IS','LO','SL','SS','UL','US'}
                    % these are int/long values
                    %fprintf('%s: Int/Long VM %d, VR %s, data = %s, val = %d\n',tag,vm,vr,tmp_data,str2num(tmp_data));
                    if (z == 1), val_data = zeros(vm,1); end
                    if (size(tmp_data,2) > 0), val_data(z) = str2double(tmp_data); end
                case {'DS','FL','FD'} 
                    % these are floating point values
                    %fprintf('%s: Float/double VM %d, VR %s, data = %s, val = %.8f\n',tag,vm,vr,tmp_data,str2num(tmp_data));
                    if (z == 1), val_data = zeros(vm,1); end
                    if (size(tmp_data,2) > 0), val_data(z) = str2double(tmp_data); end
                otherwise % just assume string
                    %error('Unknown VR = %s found!\n',vr);
                    %fprintf('Unknown VR %s, data = %s\n',vr,str_data);
                    if (z == 1), val_data = cell(vm,1); end
                    val_data{z} = tmp_data;
            end
        end
    else
        val_data = [];
    end

    junk_len = 16 * (NoOfItems - vm);
    if (NoOfItems < 1)
        junk_len = 4;
    else
        if ( (junk_len < 16) && (junk_len ~= 0) )
            junk_len = 16;
        end
    end

    %data_end_pos = ftell(fp);
    fp = fp + junk_len; % skip padding at the end of this field
    %data_final_pos = ftell(fp);

%     fprintf('%2d - ''%s''\tVM %d, VR %s, SyngoDT %d, NoOfItems %d, Data',y-1, tag, vm, vr, SyngoDT, NoOfItems);
%     if (size(str_data))
%         fprintf(' ''%s''', str_data);
%     end
%     fprintf('\n');
%     fprintf('data_len: %d  pad_len: %d  total_len: %d\n',data_end_pos-data_start_pos,data_final_pos-data_end_pos,data_final_pos-data_start_pos);
%     fprintf('VR %s\n',vr);

    hdr.(tag) = val_data;
end
