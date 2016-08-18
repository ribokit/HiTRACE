function [ areas, prof_fit, deviation ] = do_the_fit_fast( d_align, xsel, const_width, PLOT_STUFF );
% DO_THE_FIT_FAST: Fits electrophoretic traces to sums of Gaussians -- fast due to no optimization of peak positions.
%
%  [ areas, prof_fit, deviation ] = do_the_fit_fast( d_align, xsel, const_width, PLOT_STUFF );
%
%  d_align = input matrix of traces
%  xsel = band locations
%  const_width = width of gaussians [options; default = (1/4) * mean band spacing, which appears appropriate for ABI capillaries]
%  PLOT_STUFF = parameter used by GUI to suppress graphical output. [default 1]
%
% This is now a wrapper around FIT_TO_GAUSSIANS
%
% (C) R. Das 2008-2010


 [ areas, prof_fit, deviation ] = fit_to_gaussians( d_align, xsel, const_width, PLOT_STUFF );