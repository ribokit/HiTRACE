function clr = get_residue_color(seq_char)

%
% clr = GET_RESIDUE_COLOR(seq_char);
%
% Returns color code for a nucleotide. Default is magenta; blue for A,
%  green for C, orange for U/T, red for G.
%
% by T47, Rhiju Das, May 2013.
%

if ~exist('seq_char','var'); return; end;
clr = [ 1 0 1]; % magenta for unknown.

switch upper(seq_char)
    case {'A'}
        clr = [0 0 1];
    case {'C'}
        clr = [0 0.5 0];
    case {'U','T'}
        clr = [1 0.5 0];
    case {'G'}
        clr = [1 0 0];
end