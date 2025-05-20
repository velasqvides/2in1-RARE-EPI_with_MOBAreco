function [figFolder] = standardFigureDefaults()
% figure defaults

set(0,'DefaultAxesFontSize', 30)
set(0,'defaultLineMarkerSize', 9)
set(0,'defaultLineLineWidth', 2)
set(0,'defaultAxesLineWidth', 2)

%set(0,'defaultAxesFontWeight','bold')
% set(0,'defaultTextInterpreter','latex')
% set(0, 'defaultAxesTickLabelInterpreter','latex')
% set(0, 'defaultLegendInterpreter','latex');

% if ~ispc
% 
%     set(0,'defaultTextInterpreter','tex')
%     set(0, 'defaultAxesTickLabelInterpreter','tex')
%     set(0, 'defaultLegendInterpreter','tex');
% 
% end
% 
% %myFont = 'Helvetica';
% %myFont = 'DejaVu Sans';
% myFont = 'Arial';
% 
% set(0, 'defaultTextFontName', myFont)
% set(0, 'defaultAxesFontName', myFont)

figFolder = pwd;

end

