function [] = linkfig(hfigs, axstr)

    for ii=1:length(hfigs)
        figure(hfigs(ii));
        ax(ii) = gca;
    end

    % do colorbars too.
    if findstr(axstr, 'c')
        axstr(axstr == 'c') = [];
        clim = [min(min(ax.CLim)) ...
                max(max(ax.CLim))];
        for ii=1:length(ax)
            ax(ii).CLim = clim;
        end
    end

    linkaxes(ax, axstr);
end