function histo_Z_compare (Z1, Z2, counts)
% HISTO_Z_COMPARE (Z1, Z2, [counts])
%
% Plots Z score martix to histograms. Default axis [-10 10 0 100].
%
% =Input=
%   Z1, Z2          Z score matrices input, format in double.
%   [counts]        Bins in histogram. Format in double, default 2000.
%
% by T47, Mar 2013.
%

if nargin == 0; help( mfilename ); return; end;

if ~exist('counts','var'); counts=2000; end;

Zsub = Z1 - Z2;
Zplot = reshape(Zsub, 1, numel(Zsub));

figure;
hist(Zplot, counts);
axis([-10 10 0 100]);
make_lines(-0.5,'r');
