function [hash] = githash(fname, gitdir)

    if ~exist('fname', 'var')
        fname = [];
    end

    if ~exist('gitdir', 'var')
        gitdir = '';
    else
        gitdir = ['--git-dir=' gitdir];
    end

    [~, hash] = system(['TERM=xterm-256color git ' gitdir ...
                        ' log -n 1 --pretty=format:''%H'' ' fname ';']);
    % remove bash escape characters
    hash = hash(9:48)