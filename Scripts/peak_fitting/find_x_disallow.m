function [x_disallow,x_allow] = find_x_disallow( d );
%
% [x_disallow,x_allow] = find_x_disallow( d );
%


SATURATION_CUTOFF = 0.98;

x_disallow = [];
max_profile = max( d );
saturation_points = find( d >= SATURATION_CUTOFF * ...
			  max_profile );
x_disallow = union( x_disallow, saturation_points );

x_allow = setdiff( 1:length(d), x_disallow );

