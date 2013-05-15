function d_out = auto_flip_monotone(d_in, val_mnt, var_name)

%
% d_out = AUTO_FLIP_MONOTONE(d_in, val_monotone, var_name);
%
% Returns auto_flipped d array according to desired monotonicity.
%
%
% Input
% =====
%   d_in            Required        Provides the 1D array.
%   val_monotone    Optional        Provides the desired monotonicity of
%                                    d_out. Default is 1. 1 equals strictly
%                                    increasing, -1 equals strictly
%                                    decreasing.
%   var_name        Optional        Provides the name for d_in in printing
%                                    text. Default is 'd'.
%
% Output
% ======
%   d_out                           Gives the 1D array with desired
%                                    monotonicity.
%
%
% by T47, May 2013.
%

if nargin == 0; help( mfilename ); return; end;

if ~exist('val_mnt','var') || isempty(val_mnt); val_mnt = 1; end;
if (val_mnt ~= 1) && (val_mnt ~= -1); val_mnt = 1; end;
if ~exist('var_name','var') || isempty(var_name); var_name = 'd'; end;

[num_flag, str_flag] = check_monotone(d_in);
str_flag = lower(str_flag);
d_out = d_in;

fprintf(['Input ', var_name, ' (1 x ',num2str(length(d_in)),') is ', str_flag]);

if num_flag == -2 * val_mnt;
    fprintf(', FLIPPED for use.\n');
    d_out = fliplr(d_in);
elseif num_flag == 2 * val_mnt;
    fprintf(', unchanged for use.\n');
else
    fprintf(', please check.\n');
    fprintf('** Strict monotonicity required! **\n');
end;
