%      Copyright (C) 2009-2013 Biomedizinische NMR Forschungs GmbH
%                     http://www.biomednmr.mpg.de
%            Author: Tilman Johannes Sumpf <tsumpf@gwdg.de>
%
%       Distributed under the Boost Software License, Version 1.0.
%          (See accompanying file LICENSE_1_0.txt or copy at
%                 http://www.boost.org/LICENSE_1_0.txt)

classdef asSelectionClass < handle
    
    properties (GetAccess = private, SetAccess = private)
        ph        = 0;     % parent panel handle
        fh        = 0;     % parent figure handle
        dims      = [];    % dimensions of vChgArray
        
        data      = [];    % asDataClass object
        
        cmh       = struct;  % context menu handle
        
        updFig          = '';   % figure update callback
        apply2allCb     = [];   % send to all relatives callback
        apply2buttonH   = [];
        
        enabled = true ;
                
        selVcoNr= 1;     % selected vco number
        
        returnAs = 'string';   % can be either string or number
        
        defaultPlotDim = 3;
        
        initStrings = {};
        
        sendIcon = [];
    end
    
    properties (GetAccess = public, SetAccess = private)
        sendToggleState = false;
        vcos    = {};    % value changer objects
    end
    
    methods (Access = public)
        function obj = asSelectionClass(parentPanelHandle, dims, varargin)
            
            obj.ph     = parentPanelHandle;
            obj.fh     = get(obj.ph,'Parent');   % parent figure handle
            obj.dims   = dims;
            
            % create an init string
            obj.initStrings = cell(length(dims),1);
            for i = 1 : length(dims)
                obj.initStrings{i} = '1';
            end
            
            
            % evaluate varagin
            offsets = [];
            if nargin > 2
                for i=1:floor(length(varargin)/2)
                    option=varargin{i*2-1};
                    option_value=varargin{i*2};
                    switch lower(option)
                        case 'figureupdatecallback'
                            obj.updFig = option_value;
                        case 'apply2allcb'
                            obj.apply2allCb = option_value;
                        case 'initstrings'
                            obj.initStrings = option_value;
                        case 'dataobject'
                            obj.data = option_value;
                        case 'sendicon'
                            obj.sendIcon = option_value;
                        case 'offsets'
                            offsets = option_value;                            
                        otherwise
                            warning('asSelectionClass:unknownOption',...
                                'unknown option [%s]!\n',option);
                    end;
                end;
            end
            
                                    
            % create context menu
            obj.cmh.base = uicontextmenu;
            uimenu(obj.cmh.base,'Label','Enter selection string (s)'   ,...
                'callback',@(src,evnt)obj.openSetValueDlg(),...
                'Separator','on');
            
            
            % initialize array of valueChanger objects
            obj.initValueChangerArray();
            
            % consider the offset
            if ~isempty(offsets)
                obj.setOffsets(offsets, false);
            end
            
            
        end
        
        function reInit(obj, newDimensions, newSelection)
            for i = 1 : length(obj.vcos)
                obj.vcos{i}.delete(true);
            end
            delete(obj.apply2buttonH);
            
            obj.dims        = newDimensions;
            obj.initStrings = newSelection;
            
            set(obj.fh,'HandleVisibility','on');
            obj.initValueChangerArray();
            set(obj.fh,'HandleVisibility','off');
        end
        
        function tags = getColonDims(obj)
            tags = [0,0];
            for i = 1 : length(obj.vcos)
                tag = obj.vcos{i}.getColonDimTag;
                if tag ~= 0
                    tags(tag) = i;
                end
            end
        end
        
        function setColonDim1(obj, id)
            % find current colon dim 1
            if obj.enabled
                for i = 1 : length(obj.vcos)
                    tag = obj.vcos{i}.getColonDimTag;
                    if tag == 1
                        currColonDim = obj.vcos{i}.getId;
                        if currColonDim ~= id
                            % deactivate old colonDim1
                            obj.vcos{i}.setColonDimTag(0, true);
                            return;  % (don't need to check the remaining vcos)
                        end
                    end
                end
            end
        end
        
        function setColonDim2(obj, id)
            % find current colon dim 2
            if obj.enabled
                for i = 1 : length(obj.vcos)
                    tag = obj.vcos{i}.getColonDimTag;
                    if tag == 2
                        currColonDim = obj.vcos{i}.getId;
                        if currColonDim ~= id
                            % deactivate old colonDim1
                            obj.vcos{i}.setColonDimTag(0, true);
                            return;  % (don't need to check the remaining vcos)
                        end
                    end
                end
            end
        end
                
        function setPlotDim(obj, dim)
            if obj.enabled
                for i = 1 : length(obj.vcos)
                    currColonDim = obj.vcos{i}.getId;
                    if currColonDim == dim
                        obj.vcos{i}.setPlotDimTag(true);
                    else
                        obj.vcos{i}.setPlotDimTag(false);
                    end
                end
            end
        end
        
        function dim = getPlotDim(obj)
            dim = [];
            if obj.enabled
                for i = 1 : length(obj.vcos)
                    if obj.vcos{i}.getPlotDimTag
                        dim = obj.vcos{i}.getId;
                        break
                    end
                end
            end
        end
        
        function selectVco(obj,selVcoNr)
            % if no selVcoNr is given, use the last selected number from
            % the object's properties
            if obj.enabled
                if nargin == 2
                    obj.selVcoNr = selVcoNr;
                end
                
                % select the vco
                obj.vcos{obj.selVcoNr}.select;
            end
        end
        
        function ph = getPanelH(obj)
            ph = obj.vcos{obj.selVcoNr}.getPanelH;
        end
        
        function addr = getValueAsCell(obj, retNumAsStr, includeOffset)
            if nargin < 3
                includeOffset = true;
            end
            if nargin < 2
                retNumAsStr = true;
            end
            noVcos = length(obj.dims);
            addr = cell(noVcos, 1);
            for i = 1 : noVcos
                str = obj.vcos{i}.getStr(includeOffset);
                
                if retNumAsStr
                    addr{i} = str;
                else
                    [num,sflag] = str2num(str); %#ok<ST2NM> str is not necessarily a scalar
                    if sflag == 0
                        addr{i} = str;
                    else
                        addr{i} = num;
                    end
                end
            end
        end
        
        function addr = getValue(obj, includeOffset)
            if nargin < 2
                includeOffset = true;
            end
            noVcos = length(obj.dims);
            
            if strcmp(obj.returnAs, 'number')
                addr   = zeros(1,noVcos);
                
                for i = 1 : noVcos
                    addr(1,i) = obj.vcos{i}.getStr(includeOffset);
                end
            else
                if strcmp(obj.returnAs, 'string')
                    addr   = obj.vcos{1}.getStr(includeOffset);
                    if noVcos > 1
                        for i = 2 : noVcos
                            addr = strcat(addr,',',obj.vcos{i}.getStr(includeOffset));
                        end
                    end
                end
            end
        end
        
        function offs = getOffsets(obj)
            noVcos = length(obj.dims);
            offs = zeros(noVcos,1);
            for i = 1 : noVcos
                offs(i) = obj.vcos{i}.getOffset();
            end            
        end

        function setOffsets(obj, offs, runCb)
            if nargin < 3 || isempty(runCb)
                runCb = true;
            end
                
            noVcos = length(obj.dims);
            if length(offs(:)) ~= noVcos
                fprintf('Length of the offset vector has to be equal to the number of value changer (i.e.: %d)',noVcos);
                return;
            end
            for i = 1 : noVcos
                obj.vcos{i}.setOffset(offs(i), runCb);
            end            
        end
        
        function setColonDimDivisor(obj,divi)
            
            % create divisor string
            if divi == 1
                diviString = ':';
            else
                diviString = sprintf('/%d',divi);
            end
            
            % get current colon dimensions
            colDims = obj.getColonDims;
            
            % get selection string
            sel = obj.getValueAsCell(false);
            
            % if the data is at least 2D and no colon dim is given, use the
            % two first dimensions by default
            if length(sel) > 1 || sum(colDims) == 0
                colDims = [1,2];
            end
            
            % sets selection string in colon dimensions to '/2'
            for i = 1 : length(colDims)
                sel{colDims(i)} = diviString;
            end
            
            % update selection string
            obj.setValue(sel);
        end
                
        function setValue(obj,value, suppressWarning, suppressApplyToAll, suppressCallback)
            if obj.enabled
                if nargin < 5
                    suppressCallback = false;
                    if nargin < 4
                        suppressApplyToAll = false;
                        if nargin < 3
                            suppressWarning = false;
                        end
                    end
                end
                if ~iscell(value) && ~isnumeric(value)
                    value = textscan(value,'%s','Delimiter',',');
                    value = value{1};
                end
                
                lv = length(value);
                lvco = length(obj.vcos);
                l = min(lv,lvco);
                
                for i = 1 : l
                    obj.vcos{i}.setStr(value{i});
                end
                if obj.sendToggleState && ~suppressApplyToAll
                    obj.apply2allCb('selection.setValue',false, value, true, true);
                end
                if ~suppressCallback
                    obj.updFig();
                end
            end
        end
        
        function openSetValueDlg(obj)
            if obj.enabled
                prevValue = obj.getValue;
                newValue = mydlg('Enter selection string','Selection input dlg',prevValue);
                if ~isempty(newValue)
                    obj.setValue(newValue);
                end
            end
        end
        
        function increaseCurrentVc(obj)
            if obj.enabled
                obj.vcos{obj.selVcoNr}.up;
            end
        end
        
        function tag = getCurrentVcColonDimTag(obj)
            tag = obj.vcos{obj.selVcoNr}.getColonDimTag();
        end
        
        function decreaseCurrentVc(obj)
            if obj.enabled
                obj.vcos{obj.selVcoNr}.down;
            end
        end
        
        function dims = getDimensions(obj)
            dims = obj.dims;
        end
        
        function setCurrentVcValue(obj, value)
            if obj.enabled
                obj.vcos{obj.selVcoNr}.setStr(value);
                obj.updFig();
            end
        end
        
        function value = getCurrentVcValue(obj, includeOffset)
            if nargin < 2
                includeOffset = true;
            end
            if obj.enabled
                value = obj.vcos{obj.selVcoNr}.getStr(includeOffset);
            end
        end
        
        function id = getCurrentVcId(obj)
            id = obj.selVcoNr;
        end
        
        function vc = getCurrentVc(obj)
            vc = obj.vcos{obj.selVcoNr};
        end
        
        function selectNeighbour(obj, inc, refVcoId)
            if obj.enabled
                if nargin < 3
                    refVcoId = obj.selVcoNr;
                end
                
                obj.selVcoNr = refVcoId + inc;
                
                noVcos = length(obj.dims);
                
                % check if the selVco-number exists
                if obj.selVcoNr < 1
                    obj.selVcoNr = noVcos;
                else
                    if obj.selVcoNr > noVcos
                        obj.selVcoNr = 1;
                    end
                end
                
                % select the vco
                obj.vcos{obj.selVcoNr}.select;
            end
        end
        
        function toggleSend(obj, bool)
            if nargin > 1
                set(obj.apply2buttonH,'value',bool)
            end
            switch get(obj.apply2buttonH,'value')
                case 1
                    obj.sendToggleState = true;
                    obj.send()
                case 0
                    obj.sendToggleState = false;
            end
            
        end
        
        function send(obj, includeOffset)
            if nargin < 2
                includeOffset = false;
            end
            value = obj.getValue(includeOffset);
            obj.apply2allCb('selection.setValue',false, value, true, true);
        end
        
        function enable(obj, state)
            obj.enabled = state;
            set(obj.apply2buttonH,'Enable',arrShow.boolToOnOff(state));
            for i = 1 : length(obj.vcos)
                obj.vcos{i}.enable(state);
            end
        end

% defining this delete functions didn't seem to bee a good idea. Apparently
% it sometimes happens, that an orphaned object is deleted which then might
% destroy panel handles which are still in use...
%         function delete(obj)
%             if ishandle(obj.apply2buttonH)
%                 delete(obj.apply2buttonH);
%             end
%             for i = 1 : length(obj.vcos)
%                 obj.vcos{i}.delete;
%             end
%             if ishandle(obj.cmh.base)
%                 delete(obj.cmh.base);
%             end
%         end
        
    end
    
    methods (Access = private)
        
        function initValueChangerArray(obj)
            
            % preset width and height of the valueChanger
            w = .9;
            h = 2.25;
            %             border = .13;
            
            % position of the first valueChanger
            pos  = [0, 0, w, h];
            
            % number of neccessary valueChanger objects
            noVcos     = length(obj.dims);
            
            % initialize array of valueChanger
            obj.vcos    = cell(noVcos,1);
            
            % create the valueChanger objects
            colonDimCount = 0;
            for i = 1 : length(obj.vcos)
                pos(1) = (i-1) * w;
                
                colonDimTag = 0;
                if strcmp(obj.initStrings{i},':') == 1
                    colonDimCount = colonDimCount + 1;
                    if colonDimCount < 3
                        colonDimTag = colonDimCount;
                    end
                end
                
                obj.vcos{i} = asValueChangerClass(obj.ph, ...
                    'Position',pos,...
                    'max', obj.dims(i),...
                    'id',i,...
                    'initString',obj.initStrings{i},...
                    'callback',@(vcoId)obj.callback(vcoId),...
                    'KeyPressFcn',@(src, evnt)keyPressFcn(obj,src, evnt,i),...
                    'dataObject', obj.data,...
                    'colonDim1Callback',@obj.setColonDim1,...
                    'colonDim2Callback',@obj.setColonDim2,...
                    'plotdimcallback', @obj.setPlotDim,...
                    'colonDimTag', colonDimTag,...
                    'contextmenu',obj.cmh.base);
                if i == obj.defaultPlotDim
                    obj.vcos{i}.setPlotDimTag(true);
                end
            end
            
            % create send button
            pos(1) =  length(obj.vcos) * w; % left
            pos(2) =  1.73; % bottom
            pos(3:4) = [0.4264,0.4822];
            %             pos(3:4) = [0.5,0.5];
            %             pos(1) = 8.55;
            obj.apply2buttonH = uicontrol('Style','togglebutton',...
                'Parent',obj.ph,...
                'Units','centimeters',...
                'Position',pos,...
                'SelectionHighlight','off',...
                'tooltip','Send selection to relatives',...
                'Callback',@(src,evnt)obj.toggleSend(),...
                'CData',obj.sendIcon);
            
            % select first VCO by default
            if obj.selVcoNr > noVcos
                obj.selVcoNr = 1;
            end
        end
        
        function callback(obj, vcoId)
            obj.selVcoNr = vcoId;
            if obj.vcos{vcoId}.getColonDimTag && ~strcmp(obj.vcos{vcoId}.getStr,obj.vcos{vcoId}.colonStr)
                obj.vcos{vcoId}.setColonDimTag(0,true);
            end
            if obj.sendToggleState
                value = obj.getValue(false);
                obj.apply2allCb('selection.setValue',false, value, true, true);
            end
            obj.updFig();
        end
        
        function keyPressFcn(obj,~, evnt, srcVcoNr)
            % determines if left- or right arrow is pressed;
            % selects left or right neighbour of the active valueChanger
            % accordingly
            
            
            % determine the pressed key
            switch evnt.Key
                case 'leftarrow'
                    obj.selectNeighbour(-1, srcVcoNr);
                case 'rightarrow'
                    obj.selectNeighbour(+1, srcVcoNr);
            end
            
        end
    end
    
    
end
