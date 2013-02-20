function [d_out,x_warp_all, anchor_nodes] = align_by_DP_NEW( d_all, align_blocks_in, penalizeStretchFactor, slack, maxShift, windowSize,  PLOT_STUFF );
% ALIGN_BY_DP: refine alignment by non-linear warping, optimizing correlation by dynamic programming
%
%  [d_out, x_warp_all, anchor_nodes] = align_by_DP( d_all, align_blocks_in, penalizeStretchFactor, slack, maxShift, windowSize,  PLOT_STUFF );
%
% Inputs:
%  d_all = matrix with traces to be aligned
%  align_blocks_in = subsets of traces for serial alignments specified as a cell of integer vectors. Example: 
%                       specifying { [1:4], [7 5 6 8] } will first align traces 2,3, and 4 to trace 1, and
%                         then trace 5, 6, and 8 to 7.  [Default is align all to 1]
%  penalizeStretchFactor = higher values will produce local alignments that are more globally linear. [default: 10.0]
%  slack = number of pixels by which each window can expand/contract [default: 10]
%  maxShift = maximum displacement of each window boundary from its starting position [default: 100]
%  windowSize = size of window in time pixels [default: 100]
%  PLOT_STUFF = show data before and after alignments [default 1 (True)]
%  
% Outputs:
%  d_out        = matrix with aligned traces
%  x_warp_all   = [advanced] matrix describing the local realignments
%  anchor_nodes = [advanced] window boundaries. 
%
% (C) R. Das, 2009-2011, 2013
%

%if no 'block's are specified, align the whole thing to column 1
if ~exist( 'align_blocks_in' ) | length( align_blocks_in) == 0;     align_blocks_in = { [1:size(d_all,2) ] }; end
if ~iscell( align_blocks_in ); % might be a single refcol
  if length( align_blocks_in) == 1
    refcol = align_blocks_in;
    align_blocks_in = { [refcol, 1:(refcol-1), (refcol+1):size(d_all,2) ] };
  else
    align_blocks_in = { align_blocks_in };
  end
end

if ~exist( 'penalizeStretchFactor' ); penalizeStretchFactor = 10.0; end;
if ~exist( 'slack' ); slack = 10; end;
if ~exist( 'maxShift' ); maxShift = 150; end;
if ~exist( 'windowSize' ); windowSize = 100; end;
if ~exist( 'PLOT_STUFF' ); PLOT_STUFF = 1; end;

% need to double-check that the align_blocks are acceptable...
nlanes = size( d_all, 2 );
for j = 1:length( align_blocks_in )
  align_blocks_in{j} = min(max(align_blocks_in{j},1),nlanes );
end

d_out = d_all;

for j = 1:length( align_blocks_in )
  [d_align, x_warp_all, DP, choice, anchor_nodes ] = align_by_DP_block( d_out(:, align_blocks_in{j}), 1, penalizeStretchFactor, slack, maxShift, windowSize, PLOT_STUFF );
  d_out(:, align_blocks_in{j} )  = d_align;
end

if PLOT_STUFF
  colormap( 1 - gray(100 ) );
  subplot(1,2,1);
  scalefactor = 40 / mean(mean(d_all));
  image( scalefactor * d_all )
  %make_lines( [0:size(d_out,2)], 'k', 0.25 );
  
  subplot(1,2,2);
  image( scalefactor * d_out )
  %make_lines( [0:size(d_out,2)], 'k', 0.25 );
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [d_out,x_warp_all,DP,choice, anchor_nodes] = align_by_DP_block( d_all, refcol, penalizeStretchFactor, slack, maxShift, windowSize, PLOT_STUFF );

DP = [];
choice = [];

%d_ref = d_ref - mean( d_ref );
%d = d - mean( d );
which_lanes =  1:size( d_all, 2 );

num_pixels = size( d_all, 1 );

x_warp_all = zeros(  num_pixels, size( d_all, 2 ) );

% keep track of alignment 'anchor' nodes.
nodes = get_windows( windowSize, d_all(:,1) );
anchor_nodes = zeros( length( nodes )+1, size( d_all, 2) );

% parallelization! yea!
if exist( 'matlabpool' )  
  if matlabpool( 'size' ) == 0 ;   res = findResource; matlabpool( res.ClusterSize ); end
  parfor i = which_lanes;
    [x_warp_all(:,i), anchor_nodes(:,i)]    = align_by_DP_inner_loop( i, refcol, d_all, penalizeStretchFactor, slack, maxShift, windowSize, num_pixels );
end
else
  for i = which_lanes;
    [x_warp_all(:,i), anchor_nodes(:,i)]    = align_by_DP_inner_loop( i, refcol, d_all, penalizeStretchFactor, slack, maxShift, windowSize, num_pixels );
  end
end
  
d_out = apply_warp( d_all, x_warp_all, PLOT_STUFF );
   
return;  
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
function  [x_warp, anchor_nodes]  = align_by_DP_inner_loop( i, refcol, d_all, penalizeStretchFactor, slack, maxShift, windowSize, num_pixels );
  
if ( i == refcol )
  x_warp = [1:num_pixels];
  [start_nodes, final_nodes] = get_windows( windowSize, d_all(:,i) );
  anchor_nodes = [start_nodes, final_nodes(end) ];
else
  fprintf( '-- Aligning %d of %d --\n',i,size(d_all,2));
  tic
  [d_warp, x_warp, DP, choice, anchor_nodes] = refine_by_warping( (max(d_all(:,refcol),0)).^0.2, (max(d_all(:,i),0)).^0.2, penalizeStretchFactor, slack, maxShift, windowSize );    
  toc
end

x_warp = x_warp';



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ d_warp, x_warp, DP_matrix, choice, anchor_nodes ] = refine_by_warping( d_ref_input, d_ali, penalizeStretchFactor, slack, max_shift, window_size );

% window_size = 100; % in pixels
% slack = 10; % Stretch/contract each window by this many pixels
% max_shift = 100; % how far can this window drift? 

[start_nodes_ali, final_nodes_ali, NUM_WINDOWS, num_pixels] = get_windows( window_size, d_ref_input );

DP_matrix = zeros( num_pixels, NUM_WINDOWS );

% Pad reference data with zeros.
pad_size = 100;
d_ref = zeros( pad_size + num_pixels + pad_size, 1 );
relevant_middle_pixels = pad_size + [1:num_pixels];
d_ref( relevant_middle_pixels ) = d_ref_input;
%d_ref = d_ref_input;
num_pixels_pad  = length( d_ref );

norm2_ref = sum( d_ref.*d_ref );

parts_that_count = d_ref * 0;
parts_that_count( relevant_middle_pixels ) = 1;

% Precompute all necessary correlation values.
tic
start_nodes_ref_all = {};
final_nodes_ref_all = {};
conv_matrix_all = {};

for n = 1: NUM_WINDOWS

  fprintf( 'Correlation calculation for window: %d of %d\n',n, NUM_WINDOWS);

  start_node = (n-1)*window_size + 1;
  final_node = min( n*window_size, num_pixels );

  start_node = start_nodes_ali( n );
  final_node = final_nodes_ali( n );
  
  % Window of profile to be aligned.
  window_ali = [start_node : final_node];
  window_length = final_node - start_node + 1;

  min_window_length = max( window_length - slack, 2);
  max_window_length = window_length + slack;

  d_ali_sub  = d_ali( window_ali );

  % Carve out a window of the reference alignment that hopefully
  % encapsulates a segment corresponding to this piece.
  start_node_ref = max( start_node - max_shift, 1 );
  %final_node_ref = min( start_node + max_window_length + max_shift, num_pixels );
  final_node_ref = min( final_node + max_shift, num_pixels );

  num_pixels_ref_window = final_node_ref - start_node_ref + 1;

  % this is just a 'bite' of the reference profile that could potentially match 
  % the window of the profile that we are aligning. Note that it is zero padded.
  d_ref_window = zeros( pad_size + num_pixels_ref_window + pad_size );
  relevant_middle_pixels = pad_size + [1:num_pixels_ref_window];
  d_ref_window( relevant_middle_pixels ) = ...
      d_ref( pad_size +  [start_node_ref:final_node_ref ]  );
  
  parts_that_count_window = d_ref_window * 0;
  parts_that_count_window( relevant_middle_pixels ) = 1;
  
    
  %x = [1:num_pixels_pad]';
  x = [1:length(d_ref_window)]';
  block_sizes = [min_window_length:1:max_window_length];
  scales = window_length./block_sizes;
  
  d2x =interp1(1:window_length, d_ali_sub, ...
  	    x*scales, 'linear', 0.0);
  
  % should also scale intensities -- otherwise stretching windows will give a bonus?
  %d2x = d2x * (scales') ;
  
  % Do convolution calculation. Need zero padding!
  n_fft_row = length( d_ref_window );%num_pixels_pad;
  n_fft_col = length( scales );

  % Following approximates:  -1 * [profile1 - profile2 ] ^2 
  conv_matrix = 2 * real(ifft2(fftn(d_ref_window,[n_fft_row n_fft_col]) ...
			       .* fftn(d2x(end:-1:1,:), [n_fft_row n_fft_col])));

  norm2_over_window = real(ifft2(fftn( parts_that_count_window, [n_fft_row n_fft_col]) ...
				 .* fftn( d2x(end:-1:1,:).^2, [n_fft_row n_fft_col])));
  
  conv_matrix = conv_matrix - norm2_over_window;% - norm2_ref;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %  Start nodes and final nodes -- note that this retains padding offset!
  start_nodes_ref = start_node_ref - 1 + [1:length(d_ref_window)];
  
  % Disallow starting points that are outside the window (into the padding).
  % Also, disallow shifts greater than max_shift.
  goodpoints = find( start_nodes_ref > start_node_ref & ...
		     abs( start_nodes_ref - (start_node + pad_size) ) <= max_shift );

  start_nodes_ref_all{ n } = start_nodes_ref( goodpoints );

  conv_matrix_all{ n } = conv_matrix( goodpoints, :);

  [xgrid,ygrid] = meshgrid( block_sizes, start_nodes_ref( goodpoints  ));
  final_nodes_ref = xgrid + ygrid - 1;
  % Disallow ending points that are outside the window (into the padding).

  %badpoints = find( final_nodes_ref > (num_pixels_ref_window+pad_size) );
  badpoints = find( final_nodes_ref > (final_node_ref+pad_size) );

  % little check
  %s = find( scales == 1.0 );
  %[start_node final_node start_node_ref final_node_ref ]
  %[ min( start_nodes_ref(goodpoints))  min( final_nodes_ref(s,:)) final_node_ref+pad_size]
  %length( find( final_nodes_ref( s, : ) <= (final_node_ref + pad_size) ) )
  
  final_nodes_ref( badpoints ) = NaN; % push to high number

  final_nodes_ref_all{ n } = final_nodes_ref;
  
  PLOTSTUFF = 0;
  if (PLOTSTUFF)
    clf
    subplot(2,1,1)
    plot( start_node_ref:final_node_ref,d_ref_input( start_node_ref: ...
						     final_node_ref),'b' );
    hold on
    plot( start_node:final_node,d_ali( start_node: final_node),'r' );
    
    subplot(2,1,2)
    plot(  max(conv_matrix') );
    pause;
  end
  
end
toc

nodes_ali(n+1) = final_node;

%image( conv_matrix );
%pause;

% with the conventions I chose, easier to work backward, last
% window to first.
DP_matrix  = NaN * ones( num_pixels_pad, NUM_WINDOWS);
choice          = zeros( num_pixels_pad, NUM_WINDOWS );
best_block_size = zeros( num_pixels_pad, NUM_WINDOWS );

% stretch penalty. Need to make it comparable to scale of profile convolution score
PENALIZE_STRETCH = penalizeStretchFactor * ( std( d_ali ) / slack )^2;

DP_matrix( :, NUM_WINDOWS+1) = 0.0;
tic
for n = NUM_WINDOWS:-1:1 
 
  %fprintf('Dynamic programming for window: %d of %d\n',n, NUM_WINDOWS);

  start_nodes = start_nodes_ref_all{ n };
  final_nodes = final_nodes_ref_all{ n };
  conv_matrix = conv_matrix_all{n};

  for k = 1:length( start_nodes )

    in_range_points = find( ~isnan( final_nodes(k,:) ) );
    in_DP_points = find(  ~isnan( DP_matrix( final_nodes(k,in_range_points)+1, n+1 ) ) )';
    goodpoints = in_range_points( in_DP_points );
    if length( goodpoints ) > 0
      
      % New -- prevent warping if there's not much information.
      prev_score = DP_matrix( final_nodes(k, goodpoints)+1, n+1)';
      convolution_score = conv_matrix( k, goodpoints );

      stretch_score = 0 * prev_score;
      potential_block_sizes = final_nodes(k, goodpoints) - start_nodes(k) + 1;
      if ( n < NUM_WINDOWS )
	stretch_score = PENALIZE_STRETCH * (potential_block_sizes - best_block_size( final_nodes(k, goodpoints)+1, n+1 )').^2;
      end
      
      [y, i ]  = max( prev_score + convolution_score + stretch_score );      
      DP_matrix( start_nodes(k) , n) = y;      
      choice( start_nodes(k), n) = final_nodes( k, goodpoints(i) );
      best_block_size( start_nodes(k), n ) = potential_block_sizes( i );
    end

  end
  
end
toc

%clf
%image( 0.000001 * ( max(max(DP_matrix )) - DP_matrix ) )
%pause;

% Backtrack.
[y, node] = max( DP_matrix(:,1) );

start_nodes = [];
final_nodes = [];
for n = 1:NUM_WINDOWS
  start_nodes(n) = node;  
  node = choice( node, n )+1;
  final_nodes(n) = node-1;
end


start_nodes = start_nodes - pad_size;
final_nodes = final_nodes - pad_size;

nodes = [ start_nodes final_nodes(end) ];
nodes_ali = [ start_nodes_ali final_nodes_ali(end)];
%[ nodes; nodes_ali; nodes-nodes_ali]

if ( length( unique( nodes_ali )  ) ~= length( nodes_ali  ) )
  nodes_ali
  nodes
end
x_warp = interp1( nodes_ali, nodes, [1:num_pixels], 'linear','extrap');

if ( length( unique( x_warp )  ) ~= length( x_warp  ) ) % something weird.
  x_warp;
  x_warp = 1:num_pixels;
end
d_warp = interp1( x_warp, d_ali( 1: num_pixels), 1:num_pixels, 'linear',0.0);

anchor_nodes = nodes'; %useful output

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [start_nodes_ali, final_nodes_ali, NUM_WINDOWS, num_pixels] = get_windows( window_size, d_ref_input );

num_pixels = length( d_ref_input );

NUM_WINDOWS = floor( num_pixels / window_size ) + 1;

% Make sure last window is at least 1 pixel big.
if ( num_pixels - (NUM_WINDOWS-1)*window_size < 2 )
  NUM_WINDOWS = NUM_WINDOWS - 1;
end

for n = 1: NUM_WINDOWS
  start_nodes_ali(n) = (n-1)*window_size + 1;
  final_nodes_ali(n) = min( n*window_size, num_pixels );
end
