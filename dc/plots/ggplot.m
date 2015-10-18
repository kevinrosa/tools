function [] = ggplot()

    hax = gca;
    axes(hax);
    grid on;
    hax.Color = [1 1 1]*0.93;
    hax.GridColor = [1 1 1];
    hax.GridAlpha = 0.65;
    hax.TickLength = [0 0];
    hax.LineWidth = 1;
end