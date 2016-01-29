function [] = ggplot()

    hax = gca;

    grid on;
    hax.Color = [1 1 1]*0.90;
    hax.GridColor = [1 1 1];
    hax.GridAlpha = 0.65;
    hax.TickLength = [0 0];
    hax.LineWidth = 2;

    hax.XRuler.Axle.ColorData = uint8([0 0 0 0])';
    hax.YRuler.Axle.ColorData = uint8([0 0 0 0])';
end