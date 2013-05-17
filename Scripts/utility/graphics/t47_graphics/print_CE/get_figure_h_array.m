function h_array = get_figure_h_array(xsel_cell, H)

%
% h_array = GET_FIGURE_H_ARRAY(xsel_cell, H);
%
% Returns the number of sequence positions in each page row of splitted
%  xsel_cell in an array.
%
%
% Input
% =====
%   xsel_cell       Required        Provides the 3D splitted xsel positions 
%                                    and labels.
%   H               Required        Provides the number of vertical sections.
%
% Output
% ======
%   h_array                         Gives the number of sequence positions
%                                    in each page row.
%
% by T47, May 2013.
%

if nargin == 0; help( mfilename ); return; end;

h_array = zeros(1, H);
for i = 1:H
    
    xsel_sub = zeros(1, size(xsel_cell, 2));
    for j = 1:size(xsel_cell, 2)
        if ~isempty(xsel_cell{i, j, 2})
            xsel_sub(j) = xsel_cell{i, j, 2};
        end;
    end;
    h_array(i) = length(xsel_sub(xsel_sub ~= 0));
end;