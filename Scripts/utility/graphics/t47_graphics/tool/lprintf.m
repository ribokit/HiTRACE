function new_line = lprintf( prev_line, curr_line)

%
% new_line = LPRINTF (prev_line, curr_line);
%
% Show and replace a line in command window. Last line will be erased.
%
% Input
% =====
%   prev_line               The line to be erased.
%   curr_line               The line to be shown.
%
% Output
% ======
%   new_line                The line to be shown, for next cycle.
%
% by T47, Oct 2013.
%

curr_line = sprintf(curr_line);
fprintf([repmat('\b', 1, length(prev_line)), curr_line]);
new_line = curr_line;
