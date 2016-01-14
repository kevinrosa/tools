% estimate degrees of freedom using integral timescale
%        [dof, IntegralTimeScale] = calcdof(in)
function [dof, IntegralTimeScale] = calcdof(in)

    [c,lags] = xcorr(in - mean(in), 'coef');
    [cmax,imax] = max(c);

    % calculate a bunch of timescales and take maximum
    % From Talley et al., Descriptive Physical Oceanography.
    CorrArea = cumtrapz(c(imax:end));

    [~,IntegralTimeScale] = max(CorrArea./cmax);
    dof = floor(length(in)/IntegralTimeScale);
end