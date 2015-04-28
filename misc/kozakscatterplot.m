%% rangeframe
%
%  Jeffrey Boyd
%  boydjeff@gmail.com
%  October 2013
%
% This scatterplot is based on the following paper:
% Kozak, M. (2010). Improved Scatterplot Design. IEEE Computer Graphics & 
% Applications 30(6), 3-7.
%
% Link: http://doi.ieeecomputersociety.org/10.1109/MCG.2010.112
%
% The design of this scatterplot is a hybrid of designs from Edward Tufte
% and William Cleveland.
%
%% ToDo:
%
%% Function Definition
% Inputs:
%
% * *main_axis* - handle to the axis we intend to modify 
%
% Outputs:
%
% * *main_axis* - handle to the modified axis
%
function [main_axis] = kozakscatterplot(main_axis, minMaxX, minMaxY)
minX = minMaxX(1);
maxX = minMaxX(2);
minY = minMaxY(1);
maxY = minMaxY(2);
newticks = get(gca, 'XTick');
newticks(1) = minX;
newticks(end) = maxX;
while newticks(end-1) > maxX
    newticks(end-1) = [];
end
while newticks(2) < minX
    newticks(2) = [];
end
set(gca, 'XTick', newticks);
newticks = get(gca, 'YTick');
newticks(1) = minY;
newticks(end) = maxY;
while newticks(end-1) > maxY
    newticks(end-1) = [];
end
while newticks(2) < minY
    newticks(2) = [];
end
set(gca, 'YTick', newticks);
set(gca, 'TickDir', 'out');
%% get current axis properties
dims = get(main_axis, 'Position');
    left   = dims(1);           % break these out for clarity
    bottom = dims(2);           % (I wish I knew a better way...)
    width  = dims(3);
    height = dims(4);
x_axis_lims = get(main_axis, 'XLim');
    x_axis_min = x_axis_lims(1);
    x_axis_max = x_axis_lims(2);
y_axis_lims = get(main_axis, 'YLim');
    y_axis_min = y_axis_lims(1);
    y_axis_max = y_axis_lims(2);    
    
%% draw range frame

% set tick length for middle value
l_tick = .01;

% map x_frame array to normalized figure coordinates
x_frame_mapped = (minMaxX-x_axis_min)/(x_axis_max-x_axis_min) * width + left;

% map y_frame array to normalized figure coordinates
y_frame_mapped = (minMaxY-y_axis_min)/(y_axis_max-y_axis_min) * height + bottom;

% draw horizontal axis range frame
annotation('line', [x_frame_mapped(1),x_frame_mapped(2)], [bottom, bottom]);

% draw vertical axis range frame
annotation('line', [left, left], [y_frame_mapped(1),y_frame_mapped(2)]);


%% Modify the existing axes to improve the visibility of the range frame

% modify the normal axis
box on;
mask_axis = copyobj(gca, gcf); 
set(mask_axis, ...
    'color', 'none', ...
    'box', 'on', ...
    'xcolor', [.75 .75 .75], ...      % maxe the axis lines white
    'ycolor',[.75 .75 .75]', ...
    'Xtick',[], ...         % dont print any ticks or numbers
    'Ytick',[], ...
    'xgrid', 'off', ...     % make sure there is no background grid
    'ygrid','off');
grid on;


