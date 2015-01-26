function [] = center_colorbar()

    colorbar;
    clim = caxis;
    
    a = max(abs(clim));
    caxis([-a a]);

    % if this is called, then I want a diverging colorbar with
    % white at 0.
    colormap(flipud(cbrewer('div','RdBu', 32)));