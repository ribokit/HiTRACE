function [d_out, x_realign] = align_linear( d, align_blocks_in, PLOT_STUFF );
%
% ALIGN_LINEAR:  (linear-time) alignment of a matrix of electropheretic traces to first trace
%
%   [d_out, x_realign] = align_using_ref( d, align_blocks_in, PLOT_STUFF );
%
%  d     = matrix with traces to be aligned
%  align_blocks_in = subsets of traces for serial alignments specified as a cell of integer vectors. Example: 
%                       specifying { [1:4], [7 5 6 8] } will first align traces 2,3, and 4 to trace 1, and
%                         then trace 5, 6, and 8 to 7.  [Default is align all to 1]
%  Optimizes correlation coefficient by grid search + Fast Fourier Transform
%
%
% (C) R. Das, 2013
%


%if no 'block's are specified, align the whole thing to column 1
if ~exist( 'align_blocks_in' ) | length( align_blocks_in) == 0;  align_blocks_in = { [1:size(d,2) ] }; end
if ~iscell( align_blocks_in ); 
  if length( align_blocks_in) == 1% might be a single refcol
    refcol = align_blocks_in;
    align_blocks_in = { [refcol, 1:(refcol-1), (refcol+1):size(d,2) ] };
  else
    align_blocks_in = {align_blocks_in};
  end
end

if ~exist( 'PLOT_STUFF','var') PLOT_STUFF = 1; end;

[d_out, ~, x_realign ] = align_using_ref( d, d, align_blocks_in, PLOT_STUFF );