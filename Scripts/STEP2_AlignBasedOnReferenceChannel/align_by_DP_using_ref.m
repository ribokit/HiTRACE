function [d_out,d_ref_out, x_warp_all] = align_by_DP_using_ref( d, d_ref, align_blocks_in, penalizeStretchFactor, slack, maxShift, windowSize, PLOT_STUFF , SHOW_ANCHOR_NODES);
% ALIGN_BY_DP: refine alignment by non-linear warping, optimizing correlation by dynamic programming
%
% [d_out,d_ref_out, x_warp_all] = align_by_DP_using_ref( d, d_ref, align_blocks_in );
%
% (C) R. Das, 2011, 2013
%

%if no 'block's are specified, align the whole thing to column 1
if ~exist( 'align_blocks_in' ) | length( align_blocks_in) == 0;  align_blocks_in = { [1:size(d,2) ] }; end
if ~iscell( align_blocks_in ); % might be a single refcol
  refcol = align_blocks_in;
  align_blocks_in = { [refcol, 1:(refcol-1), (refcol+1):size(d,2) ] };
end

if ~exist( 'penalizeStretchFactor' ); penalizeStretchFactor = 1.0; end;
if ~exist( 'slack' ); slack = 50; end;
if ~exist( 'maxShift' ); maxShift = 200; end;
if ~exist( 'windowSize' ); windowSize = 500; end;
if ~exist( 'PLOT_STUFF' ); PLOT_STUFF = 1; end;

if ~exist( 'PLOT_STUFF' ) PLOT_STUFF = 1; end;
if ~exist( 'SHOW_ANCHOR_NODES' ) SHOW_ANCHOR_NODES = 0; end;

if ( size( d, 1) ~= size( d_ref, 1 )  |  size( d, 2) ~= size( d_ref, 2 ) ); fprintf( 'd and d_ref are not the same size!' ); end;


d_ref_out = d_ref;
d_out = d;

for j = 1:length( align_blocks_in )
  [d_ref_out, x_warp_all, anchor_nodes ] = align_by_DP( d_ref_out,  align_blocks_in{j}, penalizeStretchFactor, slack, maxShift, windowSize, PLOT_STUFF );
  d_out(:, align_blocks_in{j}) = apply_warp( d_out(:,align_blocks_in{j}), x_warp_all, 0 );
end


if PLOT_STUFF;
  colormap( 1- gray(100));

  scalefactor = 40 / mean(mean(d));  
  subplot(2,2,1);
  image( scalefactor*d );
  %make_lines( [0:size(d_out,2)], 'k', 0.25 );
  if SHOW_ANCHOR_NODES; show_anchor_nodes( anchor_nodes ); end;
  subplot(2,2,2);
  image( scalefactor*d_out );
  %make_lines( [0:size(d_out,2)], 'k', 0.25 );
  
  scalefactor_ref = 40 / mean(mean(d_ref));
  subplot(2,2,3);
  image( scalefactor_ref * d_ref );
  %make_lines( [0:size(d_out,2)], 'k', 0.25 );
  if SHOW_ANCHOR_NODES; show_anchor_nodes( anchor_nodes ); end;
  subplot(2,2,4);
  image( scalefactor_ref * d_ref_out );
  %make_lines( [0:size(d_out,2)], 'k', 0.25 );
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function show_anchor_nodes( anchor_nodes );

for i = 1:size( anchor_nodes, 1 )
  hold on
  plot( [1:size(anchor_nodes,2)], anchor_nodes(i,:), 'r-' );
end
hold off;