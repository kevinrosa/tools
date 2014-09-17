% Draw 45 degree line. takes axis handle as optional input
%     [hline] = line45(hAxis)

function [hline] = line45(hAxis)
    if ~exist('hAxis', 'var'); hAxis = gca; end

    for ii=1:length(hAxis)
        axes(hAxis(ii));
        limx = xlim;
        hline = plot(limx, limx, '--', 'Color', [1 1 1]*0.75);
    end
end