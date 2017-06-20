function cols = split_string( input_string, delimiter );
% cols = split_string( input_string, delimiter );
%
% INPUT:
%   input_string = input string
%   delimiter    = character for delimiting (\t is OK for tabs). Default: ' '.
%
% OUTPUT:
%   cols      = cell of strings that the input string was split into
%
if nargin < 1; help( mfilename); cols = {}; return; end;
if ~exist('delimiter', 'var'); delimiter = ' '; end;
if exist( 'strsplit', 'file')
    % newer versions of MATLAB have strsplit
    cols = strsplit( input_string, delimiter ); return;
end

delimiter = sprintf(delimiter);

remain = input_string;
cols = {};
while length( remain ) > 0
  [token, remain] = strtok(remain, delimiter);
  cols = [cols, token];
end


