function histo_Z (Z, counts)
% HISTO_Z (Z, [counts])
%
% Plots Z score martix to histograms. Default axis [-10 1 0 100].
%
% =Input=
%   Z               Z score matrix input, format in double.
%   [counts]        Bins in histogram. Format in double, default 1000.
%
% by T47, Mar 2013.
%

if ~exist('counts') counts=1000;    end;
Zplot = reshape(Z, 1, numel(Z));

figure;
hist(Zplot, counts);
axis([-10 1 0 100]);