function [hcbar] = center_colorbar(hcbar)

    if ~exist('hcbar', 'var')
        hcbar = colorbar;
    end
    clim = caxis;
    
    a = max(abs(clim));
    caxis([-a a]);

    % if this is called, then I want a diverging colorbar with
    % white at 0.
    colormap(flipud(cbrewer('div','RdBu', 32)));

    % always mark 0
    hcbar.Ticks = sort(unique([hcbar.Ticks 0]));
    hcbar.TicksMode = 'auto';