function [pos_sub, txt_sub] = parse_xsel_label(xsel_cell, h_sub, ct_offset)

%
% [pos_sub, txt_sub] = PARSE_XSEL_LABEL(xsel_cell, h_sub, ct_offset);
%
% Parses xsel positions and labels from the splitted cell to arrays.
% xsel_cell could be generated from xsel_label_split.
%
%
% Input
% =====
%   xsel_cell       Required        Provides the 3D splitted xsel positions 
%                                    and labels.
%   h_sub           Required        Provides the vertical page number. Only
%                                    the subset of page row h_sub will be 
%                                    in output.
%   ct_offset       Optional        Provides the count offset. Default is
%                                    1. Pairs of xsel will be output to
%                                    pos_sub and txt_sub starting with
%                                    index ct_offset. Blanks will be added
%                                    before ct_offset.
%
% Output
% ======
%   pos_sub                         Gives all the y-axis positions of xsel.
%   txt_sub                         Gives all the y-axis labels of xsel.
%
%
% by T47, May 2013.
%

if nargin == 0; help( mfilename ); return; end;

if ~exist('ct_offset','var') || isempty(ct_offset); ct_offset = 1; end;

pos_sub = []; txt_sub = {}; 

for k = 1:size(xsel_cell, 2)
    if ~isempty(xsel_cell{h_sub, k, 2});
        pos_sub(ct_offset) = xsel_cell{h_sub, k, 2};
        txt_sub{ct_offset} = xsel_cell{h_sub, k, 1};
        ct_offset = ct_offset + 1;
    end;
end;
