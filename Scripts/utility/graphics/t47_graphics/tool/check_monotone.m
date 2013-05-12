function [num_flag, str_flag] = check_monotone(d_array)

% [num_flag, str_flag] = CHECK_MONOTONE(d_array)
%
% Checks monotonicity of input array.
% d_array must be number array, not cell.
% num_flag denotes monotonicity:
%       2       Strictly monotonic increasing;
%       1       Non-strictly monotonic increasing;
%       0       Not monotonic.
%       -1      Non-strictly monotonic decreasing;
%       -2      Strictly monotonic decreasing.
% str_flag gives the answer in string.
%
% by T47, Apr 2013
%
num_flag = NaN; str_flag = '';

array_monotone = [(all(diff(d_array) <= 0)),...     % increase
    (all(diff(d_array) < 0)),...                    % strict increase
    (all(diff(d_array) >= 0)),...                   % decrease
    (all(diff(d_array) > 0))];                      % strict decrease

if array_monotone == [0 0 1 1]; num_flag = 2; str_flag = 'Strictly monotonic increasing'; end;
if array_monotone == [0 0 1 0]; num_flag = 1; str_flag = 'Non-strictly monotonic increasing'; end;
if array_monotone == [1 1 0 0]; num_flag = -2; str_flag = 'Strictly monotonic decreasing'; end;
if array_monotone == [1 0 0 0]; num_flag = -1; str_flag = 'Non-strictly monotonic decreasing'; end;
if array_monotone == [0 0 0 0]; num_flag = 0; str_flag = 'Not monotonic'; end;
