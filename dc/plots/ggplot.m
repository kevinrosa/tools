function [] = ggplot()

    hax = gca;
    axes(hax);
    grid on;
    hax.Color = [1 1 1]*0.95;
    hax.GridColor = [1 1 1];
    hax.GridAlpha = 0.65;
end