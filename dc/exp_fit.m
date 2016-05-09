% I want to fit y = y_1 + y_0 exp(-(x-x0)/X)
function [y0, X, x0, y1] = exp_fit(x, y, plot_flag, test)

    if ~exist('test', 'var'), test = 0; end
    if ~exist('plot_flag', 'var'), plot_flag = 0; end

    if test
        test_exp_fit();
        return;
    end

    x = double(x); y = double(y);
    x = x(~isnan(y)); y = y(~isnan(y));

    if size(y) ~= size(x), y = y'; end

    initGuess(1) = max(abs(y(:)));
    initGuess(2) = mean(x(:));
    initGuess(3) = mean(y(:));
    initGuess(4) = 0;

    opts = optimset('MaxFunEvals',1e7, 'TolFun', 1e-10);
    [fit2,~,exitflag] = fminsearch(@(fit) fiterror(fit,x,y), ...
                                   initGuess,opts);

    y0 = fit2(1);
    X = fit2(2);
    y1 = fit2(3);
    x0 = fit2(4);

    if plot_flag
        figure;
        plot(x,y,'k*'); hold all
        plot(x,y1 + y0*exp(-(x-x0)/X));
    end
end

function [E] = fiterror(fit,x,y)
    % x = (T0,H,a)
    y0 = fit(1); X = fit(2); y1 = fit(3); x0 = fit(4);

    E = sum((y - y1 - y0.*exp(-(x-x0)/X)).^2);
end

function [] = test_exp_fit()
    x = [0.02:0.05:5];
    X = 2;
    y0 = 2;
    y1 = 10;
    x0 = 2;
    y = y1 + y0 * exp(-(x-x0)/X) + 0; rand(size(x)) .* y0/4;

    [Y0, XX, Y1, X0] = exp_fit(x,y,1);

    disp(['y0 = ' num2str(Y0) ' | Original = ' num2str(y0)]);
    disp(['X = ' num2str(XX) ' | Original = ' num2str(X)]);
    disp(['y1 = ' num2str(Y1) ' | Original = ' num2str(y1)]);
    disp(['x0 = ' num2str(X0) ' | Original = ' num2str(x0)]);
end