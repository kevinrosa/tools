%% Kozak Scatterplot Demo
%  
%  Jeffrey Boyd
%  boydjeff@gmamil.com
%  October 2013
%
% Requires kozakscatterplot.m
%
% This demo creates a scatterplot as described by Marcin Kozak in
% Kozak, M. (2010). Improved Scatterplot Design. IEEE Computer Graphics & 
% Applications 30(6), 3-7.
%
% Link: http://doi.ieeecomputersociety.org/10.1109/MCG.2010.112
%
% The design of this scatterplot is a hybrid of designs from Edward Tufte
% and William Cleveland.
%
figure;
%% Generate data to plot
x = exp(0 + .5*randn(1,100));                 % the log-normal is an interesting distribution
y = x + .15*randn(1,100);

minMaxX = [min(x) max(x)];
minMaxY = [min(y) max(y)];

%% plot the data
% we can use a simple plot command, as the changes are made at the axis level
% instead of at the figure level.
plot(x, y, '.');

% apply a range frame to the data we are interested in
ax1 = kozakscatterplot(gca, minMaxX, minMaxY);

% we can make changes to the axis after it is returned from the rangeframe function
xlabel(ax1, 'X data');                          %label everything
ylabel(ax1, 'Y data');
legend show;
title(ax1, 'Statistical data in a Kozak Scatterplot');