% Makes all graphs kickass. Call after making all modifications to plot
%       function [] = beautify(fontSizes)
% fontSizes = [Axis,Labels,Title] = [12 12 14] - default

%TODO: Support for multiple axes

% Deepak Cherian 05/11/2010

function [] = beautify(fontSizes)

    drawnow;

    if ~exist('fontSizes','var') || isempty(fontSizes)
        fontSizes = [16 16 18];
    end

    if strcmpi(computer,'GLNXA64')
        font_name = 'Helvetica';
    else
        font_name = 'TeXGyrePagella'; %'AvantGarde'
    end
    font_name_axis = font_name;

    % Get required handles for current figure
    hFig = evalin('caller','gcf');
    hAxis = evalin('caller','gca');

    if ~ishandle(hFig) || ~ishandle(hAxis)
     fprintf('\n Error in beautify.m. Invalid Handle');
    end

    % Get some more handles
    hXLabel = get(hAxis,'XLabel');
    hYLabel = get(hAxis,'YLabel');
    hZLabel = get(hAxis,'ZLabel');
    hTitle  = get(hAxis,'Title');

    % Aaaand..  Presto!
    set(hFig, ...
        'Color','white', ...
        'renderer', 'painters');
    set(legend, ...
        'box', 'off');

    set(hTitle, ...
        'FontSize'   , fontSizes(3) , ...
        'FontWeight' , 'normal', ...
        'FontName'   , font_name);

    set([hXLabel, hYLabel, hZLabel]  , ...
        'FontName'   , font_name, ...
        'FontWeight' , 'normal', ...
        'FontSize'   , fontSizes(2)         );

    % rotate y-label and make space for it
    hYLabel.Rotation = 0;
    pos = hYLabel.Position;
    ex = hYLabel.Extent;
    hYLabel.Position(1) = pos(1) - ex(3)/2;

    % y-label will now be outside figure space.
    % need to reduce axis width

    set(hAxis, ...
        'FontName'    , font_name_axis, ...
        'Box'         , 'on'     , ...
        'TickDir'     , 'in'     , ...
        'TickLength'  , [.005 .005] , ... % IMPROVE THIS
        'FontWeight'  , 'normal', ...
        ...%'XMinorTick'  , 'on'      , ...
        ...%'YMinorTick'  , 'on'      , ...
        ...%'ZMinorTick'  , 'on'      , ...
        'XGrid'       , 'off'      , ...
        'YGrid'       , 'off'      , ...
        'ZGrid'       , 'off'      , ...
        'XColor'      , [.3 .3 .3], ...
        'YColor'      , [.3 .3 .3], ...
        'ZColor'      , [.3 .3 .3], ...
        'FontSize'    , fontSizes(1), ...
        'LineWidth'   , 2        );

     % Line Width 2
     set(findobj('Type','line'),'LineWidth',2)

     % find contours / images and then set box on + renderer = zbuffer
     if ~isempty(findall(gcf,'type','contourgroup','visible','on')) || ...
               ~isempty(findall(gcf,'type','image','visible','on'))
         set(hFig,'renderer','zbuffer');
         setset(hAxis,'box','on');
     end

%      warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
%      jframe = get(hFig,'JavaFrame');
%      jframe.setMaximized(true);

    % Makes Fullscreen
   %set(hFig,'un','n','pos',[0,0,0.95,0.95]); figure(hFig);

%    %Remove whitespace around graph - use export_fig
%     T = get(hAxis,'tightinset');
%     set(hAxis,'position',[T(1) T(2) 1-T(1)-T(3) 1-T(2)-T(4)]);

    % WYSIWYG output apparently if print called after beautify
    %set(hFig, 'PaperPositionMode', 'auto');