% Correct overlapping tick labels
% Give it an axis and which tick index to remove
% also control number of decimal places using
% C style format string
%
%        correct_ticks(axis_name, format_string, tick_indices)
function [] = correct_ticks(ax, format, ind)

    tickstr = [upper(ax) 'Tick'];
    labelstr = [tickstr 'Label'];

    ticks = get(gca, tickstr);
    if ~ischar(ind) && ~iscell(ind) % provided with index
        ind(ind > length(ticks)) = [];
        ticks(ind) = [];
    else
        if ~iscell(ind)
            ind = cellstr(ind);
        end
        for kk=1:length(ind)
            ticks(ticks == str2double(ind{kk})) = [];
        end
    end
    set(gca, tickstr, ticks);

    if ~isempty(format)
        for ii = 1:length(ticks)
            ticklabels{ii} = sprintf(format, ticks(ii));
            if ticks(ii) == 0
                ticklabels{ii} = '0';
            end
        end
        set(gca, labelstr, ticklabels);
    end
end