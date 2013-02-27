function [ymin, ymax] = findTimeRange(d)
% findTimeRange
%
% [ymin, ymax] = findTimeRange(d)
%
% Input: 
%  d = set of traces for which we need to find edges
%
% Output:
%
%  ymin = minimum time pixel for starting edge
%  ymax = maximum time pixel for final edge
%
% by kprotoss
% find time range automatically
% modified a little by Rhiju based on some of his data sets.
% 

if nargin == 0;  help( mfilename ); return; end;

N = size(d,2);

for j = 1:N;
    y = edge(d(:,j));
    y = find( y ~= 0);

    min_list(j) = y(1);
    max_list(j) = y(end);
end

% get rid of outliers...
ymin = round( min_list( ceil( 0.50*N ) ) - 200);
ymax = round( max_list( ceil( 0.90*N ) ) + 200);

ymin = max( ymin, 1 );
ymax = min( ymax, size(d,1) );