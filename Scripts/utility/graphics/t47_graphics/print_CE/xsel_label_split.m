function [xsel_all, xsel_split] = xsel_label_split(xsel, seqpos, sequence, offset, H, h_length, h_sp)

%
% [xsel_all, xsel_split] = XSEL_LABEL_SPLIT(xsel, seqpos, sequence, offset, ...
%                                           H, h_length, h_sp);
%
% Splits xsel to 1 x H sub matrices for each figure row.
% H, h_length could be generated from d_expand_divisible.
%
%
% Input
% =====
%   xsel            Required        Provides the xsel positions.
%   seqpos          Required        Provides the sequence position numberings.
%                                    Used fo each nucleotide's numbering.
%   sequence        Required        Provides the sequence. Used for each
%                                    nucleotide's labeling.
%   offset          Required        Provides the sequence numbering offset.
%   H               Required        Provides the number of vertical sections.
%   h_length        Required        Provides the number of rows in each section.
%   h_sp            Required        Provides the vertical spacer size to add on
%                                    each side of each section.
%
% Output
% ======
%   xsel_all                        Gives all xsel positions and corresponding
%                                    labels. This is a 2D cell, with 1st column
%                                    the text string in 'G145' format, 2nd
%                                    column xsel y-axis coordinates.
%   xsel_split                      Gives xsel positions and corresponding
%                                    labels splitted to each figure row. This
%                                    is a 3D cell, with 1st dimension the
%                                    vertical page number, 2nd and 3rd same as
%                                    in xsel_all.
%
%
% by T47, May 2013.
%

if nargin == 0; help( mfilename ); return; end;

% read in band names (Y-axis)
xsel_all = cell(length(seqpos), 2);
for i = 1:length(seqpos)
    xsel_all{i, 1} = [sequence(seqpos(i) - offset), num2str(seqpos(i))];
    xsel_all{i, 2} = xsel(i);
end;

% split band position array into h*h_length array
xsel_split = cell(H, size(xsel_all, 1), 2);
for i = 1:H
    ymin = h_length * (i - 1) + 1; ymax = h_length * i;
    for j = 1:size(xsel_all, 1)
        if xsel_all{j, 2}>=ymin && xsel_all{j, 2}<=ymax;
            xsel_split{i, j, 1} = xsel_all{j, 1}; 
            xsel_split{i, j, 2} = xsel_all{j, 2} + h_sp - (ymin - 1);
        end;
    end;
end;
