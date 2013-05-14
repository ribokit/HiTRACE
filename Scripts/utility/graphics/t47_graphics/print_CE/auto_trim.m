function [d_trim, xsel_trim] = auto_trim(d_org, xsel_org, up_offset, low_offset, flg)

%
% [d_trim, xsel_trim] = AUTO_TRIM(d_org, xsel_org, up_offset, low_offset, flag);
%
% Returns trimmed d_align and shifted xsel based up_offset and low_offset. 
%
%
% Input
% =====
%   d_org           Required        Provides the original d matrix.
%   xsel_org        Required        Provides the original xsel positions.
%   up_offset       Optional        Provides the excess area of d before the first
%                                    xsel position. Default is 50.
%   low_offset      Optional        Provides the excess area of d after the last 
%                                    xsel position. Default is 50.
%   flag            Optional        Provides whether to perform auto_trim.
%                                    Default is 1 (YES). 0 to disable auto_trim 
%                                    and take original d.
%
% Output
% ======
%   d_trim                          Gives the trimmed d matrix.
%   xsel_trim                       Gives the shifted xsel positions.
%
%
% by T47, May 2013.
%

if nargin == 0; help( mfilename ); return; end;

if ~exist('up_offset','var') || isempty(up_offset); up_offset = 50; end;
if ~exist('low_offset','var') || isempty(low_offset); low_offset = 50; end;
if ~exist('flg','var') || isempty(flg); flg = 1; end;

d_upper_bound = max([round((min(xsel_org) - up_offset)), 1]);
d_lower_bound = min([round((max(xsel_org) + low_offset)), size(d_org, 1)]);

if flg == 1;
    d_trim = d_org(d_upper_bound:d_lower_bound, :);
    xsel_trim = xsel_org - d_upper_bound + 1;
    
    fprintf([', with upper offset (', num2str(up_offset), ') and lower offset (', ...
        num2str(low_offset), '), and trimmed to ', num2str(d_upper_bound), ' : ', num2str(d_lower_bound)]);
end;
