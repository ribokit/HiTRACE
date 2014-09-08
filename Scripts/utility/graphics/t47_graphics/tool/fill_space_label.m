function str_out = fill_space_label(str_in, pos)

%
% str_out = FILL_SPACE_LABEL(str_in, pos);
%
% Returns a set of strings with same length that filled up with spaces of
%  the input strings.
% e.g. {'aa', 'A'} will get {'aa', ' A'}.
%
% Input
% =====
%   str_in          Required        Provides the input string cell.
%   pos             Optional        Provides whether to add the spaces to
%                                    the end. Default is 1. 0 means add 
%                                    spaces in the head.
%
% Output
% ======
%   str_out                         Gives the output string cell with same
%                                    length of each element.
%
% by T47, May 2013.
%

if nargin == 0; help( mfilename ); return; end;

if ~exist('pos','var') || isempty(pos); pos = 1; end;
max_l = max(cellfun(@length, str_in));
str_out = cell(size(str_in));
for i = 1:length(str_in)
    if pos;
        str_out{i} = [str_in{i} blanks(max_l - length(str_in{i}))];
    else
        str_out{i} = [blanks(max_l - length(str_in{i})) str_in{i}];
    end;
end;
