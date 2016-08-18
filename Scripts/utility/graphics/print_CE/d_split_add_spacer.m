function d_spl = d_split_add_spacer(d_org, H, h_length, h_sp, W, w_length, w_sp)

%
% d_spl = D_SPLIT_ADD_SPACER(d_org, H, h_length, h_sp, W, w_length, w_sp);
%
% Splits a d matrix into H x W sub matrices with spacers added on each border.
% H, h_length, W, w_length could be generated from d_expand_divisible.
%
%
% Input
% =====
%   d_org       Required        Provides the original d matrix, expanded to
%                                be divisible by H and W.
%   H           Required        Provides the number of vertical sections.
%   h_length    Required        Provides the number of rows in each section.
%   h_sp        Optional        Provides the vertical spacer size to add on
%                                each side of each section. Default is 50.
%   W           Required        Provides the number of horizontal sections.
%   w_length    Required        Provides the number of columns in each section.
%   w_sp        Optional        Provides the horizontal spacer size to add on
%                                each side of each section. Default is 1.
%
% Output
% ======
%   d_spl                       Gives the splitted d matrices. This is a 4D
%                                matrix, with 1st dimension horzontal page
%                                number, 2nd dimension vertical page
%                                number, 3rd and 4th dimension same as
%                                d_org: e.g. d_spl(1, 2, :, :) will be the
%                                sub matrix of second page on first row.
%
%
% by T47, May 2013.
%

if nargin == 0; help( mfilename ); return; end;

if ~exist('h_sp','var') || isempty(h_sp); h_sp = 50; end;
if ~exist('w_sp','var') || isempty(w_sp); w_sp = 1; end;

d_spl = zeros(W, H, (h_length + 2 * h_sp), (w_length + 2 * w_sp)); 

for i = 1:W
    for j = 1:H
        d_spl(i, ...
              j, ...
              ((h_sp + 1):(h_sp + h_length)), ...
              ((w_sp + 1):(w_sp + w_length)) ) ...
            = d_org((h_length * (j - 1) + 1):(h_length * j),...
                    (w_length * (i - 1) + 1):(w_length * i));
    end;
end;
