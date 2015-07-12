function [] = linkfig(hfigs, axstr)

    for ii=1:length(hfigs)
        figure(hfigs(ii));
        ax(ii) = gca;
    end

    linkaxes(ax, axstr);
end