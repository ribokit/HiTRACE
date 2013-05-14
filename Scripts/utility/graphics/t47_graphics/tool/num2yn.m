function str = num2yn(num)

% str = NUM2YN(num)
%
% Returns string "YES" or "NO" (upper case) from a boolean-double input (0 or 1).
% Returns "ERR" if input is invalid.
%
% by T47, Apr 2013
%

if nargin == 0; help( mfilename ); return; end;

if num == 1; 
    str = 'YES'; 
elseif num == 0;
    str = 'NO'; 
else
    str = 'ERR';
end;
