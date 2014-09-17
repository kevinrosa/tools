% Draws horizontal or vertical line and adds tick mark
% called by linex and liney
function [handles, txthandles] = dcline(ax,x,label,color)

    hFig = evalin('caller','gcf');
    hAxis = evalin('caller','gca');
    
    if size(x,2) == 1, x = x'; end
    
    figure(hFig);
    hold on;
    
    if ~exist('label','var'), label = [];num2str(x'); end
    if ~exist('color','var') || isempty(color), color = [1 1 1]*0.7; end
    
    if length(x) ~= size(label,1)
        label = repmat(label,[length(x) 1]);
    end
    
    if ~iscell(label) && ~isempty(label)
        label = cellstr(label);
    end
    
    if ax == 'x'
        yax  = get(hAxis,'YLim');
        tickstr = 'XTick';
    else
        xax  = get(hAxis,'XLim');
        tickstr = 'YTick';
    end
    
    tick = get(hAxis,tickstr);
    
    handles = nan(size(x));
    txthandles = nan(size(x));
    
    for i=1:length(x)
        if ax == 'x'
            handles(i) = plot([x(i) x(i)],yax,'--','LineWidth',2,'Color',color);
            if ~isempty(label)
                tick = get(hAxis, 'YTick');
                txthandles(i) = text(double(x(i)),double(tick(end-1)),label{i}, ...
                    'Rotation',90,'VerticalAlignment','Bottom','FontSize',16,'Color',color);%,'FontWeight','Bold');
            end
        else
            handles(i) = plot(xax,[x(i) x(i)],'--','LineWidth',2,'Color',color);
            if ~isempty(label)
                tick = get(hAxis, 'XTick');
                txthandles(i) = text(double(tick(end-1)),double(x(i)),label{i}, ...
                    'Rotation',0,'VerticalAlignment','Bottom','FontSize',16,'Color',color);%,'FontWeight','Bold');
            end
        end
    end
    
    % add extra axis ticks
    %set(hAxis,tickstr,sort(unique(str2num(sprintf('%.2e ', [tick x])))));
