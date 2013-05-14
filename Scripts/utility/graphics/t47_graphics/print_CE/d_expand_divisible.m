function [d_exp, h_length, w_length] = d_expand_divisible(d_org, H, W)

%
% [d_exp, h_length, w_length] = D_EXPAND_DIVISIBLE(d_org, H, W);
%
% Expands d matrix to divisible sizes according to provided H and W
%  dimensions. Blank rows and/or columns will be added in the end.
%
%
% Input
% =====
%   d_org       Required        Provides the original d matrix.
%   H           Required        Provides the number of vertical sections.
%   W           Required        Provides the number of horizontal sections.
%
% Output
% ======
%   d_exp                       Gives the expanded d matrix with blank rows
%                                and/or columns added.
%   h_length                    Gives the number of rows in each section.
%   w_length                    Gives the number of columns in each section.
%
%
% by T47, May 2013.
%

if nargin == 0; help( mfilename ); return; end;

d_exp = d_org;

h_length = floor(size(d_org, 1)/H);
if h_length * H ~= size(d_org, 1); h_length = h_length + 1; end;
d_exp((size(d_org, 1) + 1):(h_length * H), :) = 0;

w_length = floor(size(d_org, 2)/W);
if w_length * W ~= size(d_org, 2); w_length = w_length + 1; end;
d_exp(:, (size(d_org, 2) + 1):(w_length * W)) = 0;