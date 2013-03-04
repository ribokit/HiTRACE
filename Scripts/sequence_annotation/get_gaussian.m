function G = get_gaussian(x, center, width )
% G = get_gaussian(x, center, width )
%
% G = exp( - ( x - center ).^2 / (2*width^2) );
%

if nargin == 0;  help( mfilename ); return; end;

G = exp( - ( x - center ).^2 / (2*width^2) );