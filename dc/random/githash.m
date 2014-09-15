function [hash] = githash(fname)

    if ~exist('fname', 'var')
        fname = [];
    end

    [~, hash] = system(['TERM=xterm-256color git log -n 1 ' ...
                        '--pretty=format:''%H'' ' fname ';']);
    % remove bash escape characters
    hash = hash(9:48)