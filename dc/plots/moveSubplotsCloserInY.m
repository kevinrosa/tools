function [] = moveSubplotsCloserInY(m,n,hax)
    Nax = length(hax);

    for ii=1:n
        dy = hax(ii).OuterPosition - hax(ii+n).OuterPosition;
        dy = dy(2);

        hax(ii).OuterPosition(2) = hax(ii).OuterPosition(2) - dy/8;
        hax(ii+n).OuterPosition(2) = hax(ii+n).OuterPosition(2) + dy/8;
    end

end