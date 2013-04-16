function [d_out, d_ref_out, x_realign] = align_using_ref( d, d_ref, align_blocks_in, PLOT_STUFF );
%
% ALIGN_USING_REF:  (linear-time) alignment of a matrix of electropheretic traces to first trace
%
%   [d_out, d_ref_out, x_realign] = align_using_ref( d, d_ref, align_blocks_in, PLOT_STUFF );
%
%  d     = matrix with traces to be aligned
%  d_ref = matrix with reference traces, which will determin alignment
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


d_ref_out = d_ref;
d_out = d;

N = size( d, 1 );
for j = 1:length( align_blocks_in )

  align_block = align_blocks_in{j};
  
  [d_ref_out(:,align_block), x_realign] = align_to_first_ver3( d_ref_out(:,align_block) );
  for i = 1:length(align_block)
    d_out(:,align_block(i))  = interp1( [1:N]', d_out(:,align_block(i)), x_realign(:,i), 'linear',0.0 ); 
  end

end


if PLOT_STUFF
  colormap( 1- gray(100));

  scalefactor = 40 / mean(mean(d));  
  subplot(2,2,1);
  image( scalefactor*d );
  subplot(2,2,2);
  image( scalefactor*d_out );
  
  scalefactor_ref = 40 / mean(mean(d_ref));
  subplot(2,2,3);
  image( scalefactor_ref * d_ref );
  subplot(2,2,4);
  image( scalefactor_ref * d_ref_out );
end

  

