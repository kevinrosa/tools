function [] = linkfig(hfigs, axstr)

    for ii=1:length(hfigs)
        figure(hfigs(ii));
        ax(ii) = gca;
    end

    % do colorbars too.
    if findstr(axstr, 'c')
        axstr(axstr == 'c') = [];
        clim = ax(1).CLim;
        for ii=1:length(ax)
            clim(1) = max([clim(1) ax(ii).CLim(1)]);
            clim(2) = min([clim(2) ax(ii).CLim(2)]);
        end
        for ii=1:length(ax)
            ax(ii).CLim = clim;
        end
    end

    linkaxes(ax, axstr);
end