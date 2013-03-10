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

min_list = sort_and_filter( min_list );
max_list = sort_and_filter( max_list );

%subplot(2,1,1);
%plot( min_list );
%subplot( 2,1,2);
%plot( max_list );
%pause;

% get rid of outliers...
N = length( min_list );
ymin = round( min_list( ceil( 0.50*N ) ) ) - 150;

N = length( max_list );
ymax = round( max_list( ceil( 0.80*N ) ) ) + 150;

ymin = max( ymin, 1 );
ymax = min( ymax, size(d,1) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function min_list = sort_and_filter( min_list );

min_list = sort( min_list );

N = length( min_list );
q1 = min_list( ceil( 0.25*N) );
q3 = min_list( ceil( 0.75*N) );

gp = find(  min_list >= (q1 - 1.5 *abs( q3-q1)) & ...
	    min_list <= (q3 + 1.5 *abs( q3-q1)) );
min_list = min_list( gp );


