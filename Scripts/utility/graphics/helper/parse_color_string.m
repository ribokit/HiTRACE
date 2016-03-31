function clr_cell = parse_color_string(clr_str)

%
% clr_cell = PARSE_COLOR_STRING(clr_str);
%
% Returns a string cell contains each color letter.
%
% by T47, May 2013.
%

clr_cell = {};
if ~exist('clr_str','var') || isempty(clr_str); return; end;

for i = 1:length(clr_str)
    clr_cell{i} = clr_str(i);
end;
