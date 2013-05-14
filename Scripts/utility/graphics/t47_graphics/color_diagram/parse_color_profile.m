function [color_scheme, d_offset, max_color, min_color] = parse_color_profile(color_profile)

%
% [color_scheme, d_offset, max_color, min_color] = PARSE_COLOR_PROFILE(color_profile) 
%
% Parses color_profile to 4 variables.
%
%
% by T47, May 2013.
%

if nargin == 0; help( mfilename ); return; end;

color_scheme = color_profile(1);
d_offset = color_profile(2);
max_color = color_profile(3);
min_color = color_profile(4);