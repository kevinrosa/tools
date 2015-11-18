% estimate degrees of freedom using integral timescale
%        [dof, IntegralTimeScale] = calcdof(in)
function [dof, IntegralTimeScale] = calcdof(in)

    [c,lags] = xcorr(in - mean(in), 'coef');
    [cmax,imax] = max(c);
    i0 = imax + find(c(imax:end) < 0, 1, 'first') - 1; % first zero crossing

    % calculate a bunch of timescales and take maximum
    % From Talley et al., Descriptive Physical Oceanography.
    for tt=imax+3:length(c)
        ff = tt - (imax+3) + 1;
        CorrArea = trapz(c(imax:tt));
        IntegralTimeScale(ff) = CorrArea/cmax;
    end

    IntegralTimeScale = max(IntegralTimeScale);
    dof = floor(length(in)/IntegralTimeScale);
end