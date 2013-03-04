function Q=quartiles(X)
%   Q=quartiles(X)
if nargin == 0;  help( mfilename ); return; end;

Nx = size(X,1);

% compute mean
mx = mean(X);

% compute the standard deviation
sigma = std(X);

% compute the median
medianx = median(X);

% STEP 1 - rank the data
y = sort(X);

% compute 25th percentile (first quartile)
Q(1) = median(y(find(y<median(y))));

% compute 50th percentile (second quartile)
Q(2) = median(y);

% compute 75th percentile (third quartile)
Q(3) = median(y(find(y>median(y))));
