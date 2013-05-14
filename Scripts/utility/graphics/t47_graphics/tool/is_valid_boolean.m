function out = is_valid_boolean (in)

% out = IS_VALID_BOOLEAN(in);
%
% Check validity of boolean (0/1) input.
%
% Input
% =====
%   in      Required        Provides a number boolean input.
%
% Output
% ======
%   out                     Provides a valid boolean output. Non-0-or-1
%                               numbers will be regard as 0. A warning will
%                               appear on console.
%
% by T47, May 2013.
%

if nargin == 0; help( mfilename ); return; end;

out = in;
if (in ~= 0) && (in ~= 1);
    out = 0;
    fprintf('WARNING: Invalid boolean flag input.\n');
end;