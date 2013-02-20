function [d_out,x_warp_all,DP,choice] = align_by_DP( d_all, align_blocks_in, slack, maxShift, windowSize,  PLOT_STUFF );
% ALIGN_BY_DP: refine alignment by non-linear warping, optimizing correlation by dynamic programming
%
%  [d_out,x_warp_all,DP,choice] = align_by_DP( d_all, slack, maxShift, windowSize, align_blocks_in, PLOT_STUFF );
%
% (C) R. Das, 2009-2011
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
if ~exist( 'PLOT_STUFF' ); PLOT_STUFF = 1; end;
if ~exist( 'slack' ); slack = 10; end;
if ~exist( 'windowSize' ); windowSize = 100; end;
if ~exist( 'maxShift' ); maxShift = 100; end;

% need to double-check that the align_blocks are acceptable...
nlanes = size( d_all, 2 );
for j = 1:length( align_blocks_in )
  align_blocks_in{j} = min(max(align_blocks_in{j},1),nlanes );
end

d_out = d_all;

for j = 1:length( align_blocks_in )
  [d_align, x_warp_all, DP, choice ] = align_by_DP_block( d_out(:, align_blocks_in{j}), 1, slack, maxShift, windowSize, PLOT_STUFF );
  d_out(:, align_blocks_in{j} )  = d_align;
end


if PLOT_STUFF
  colormap( 1 - gray(100 ) );
  subplot(1,2,1);
  scalefactor = 40 / mean(mean(d_all));
  image( scalefactor * d_all )
  
  subplot(1,2,2);
  image( scalefactor * d_out )
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [d_out,x_warp_all,DP,choice] = align_by_DP_block( d_all, refcol, slack, maxShift, windowSize, PLOT_STUFF );

DP = [];
choice = [];

%d_ref = d_ref - mean( d_ref );
%d = d - mean( d );
which_lanes =  1:size( d_all, 2 );

num_pixels = size( d_all, 1 );

x_warp_all = zeros(  num_pixels, size( d_all, 2 ) );

% parallelization! yea!
if exist( 'matlabpool' )  
  if matlabpool( 'size' ) == 0 ;   res = findResource; matlabpool( res.ClusterSize ); end
  parfor i = which_lanes;
    x_warp_all(:,i)    = align_by_DP_inner_loop( i, refcol, d_all, slack, maxShift, windowSize, num_pixels );
  end
else
  for i = which_lanes;
    x_warp_all(:,i)    = align_by_DP_inner_loop( i, refcol, d_all, slack, maxShift, windowSize, num_pixels );
  end
end
  
d_out = apply_warp( d_all, x_warp_all, PLOT_STUFF );
   
return;  
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
function   x_warp    = align_by_DP_inner_loop( i, refcol, d_all, slack, maxShift, windowSize, num_pixels );
  
if ( i == refcol )
  x_warp = [1:num_pixels];
else
  fprintf( '-- Aligning %d of %d --\n',i,size(d_all,2));
  tic
  [d_warp, x_warp, DP, choice] = refine_by_warping( (abs(d_all(:,refcol))).^0.2, (abs(d_all(:,i))).^0.2, slack, maxShift, windowSize );    
  toc
end

x_warp = x_warp';



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ d_warp, x_warp, DP_matrix, choice ] = refine_by_warping( d_ref_input, d_ali, slack, max_shift, window_size );

% window_size = 100; % in pixels
% slack = 10; % Stretch/contract each window by this many pixels
% max_shift = 100; % how far can this window drift? 

num_pixels = length( d_ref_input );
NUM_WINDOWS = floor( num_pixels / window_size ) + 1;

% Make sure last window is at least 1 pixel big.
if ( num_pixels - (NUM_WINDOWS-1)*window_size < 2 )
  NUM_WINDOWS = NUM_WINDOWS - 1;
end

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

  start_nodes_ali( n  ) = start_node;
  final_nodes_ali( n  ) = final_node;

  % Window of profile to be aligned.
  window_ali = [start_node : final_node];
  window_length = final_node - start_node + 1;

  min_window_length = max( window_length - slack, 2);
  max_window_length = window_length + slack;

  d_ali_sub  = d_ali( window_ali );

  % Carve out a window of the reference alignment that hopefully
  % encapsulates a segment corresponding to this piece.
  start_node_ref = max( start_node - max_shift, 1 );
  final_node_ref = min( start_node + max_window_length + max_shift, num_pixels );
  num_pixels_ref_window = final_node_ref - start_node_ref + 1;

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
  %  Start nodes and final nodes -- retain padding offset
  start_nodes_ref = start_node_ref - 1 + [1:length(d_ref_window)];

  % Disallow starting points that are outside the window (into the padding).
  % Also, disallow shifts greater than max_shift.
  %  goodpoints = find( start_nodes_ref > start_node_ref & ...
  %		     abs( start_nodes_ref - start_node ) <= max_shift );
  goodpoints = find( start_nodes_ref > start_node_ref & ...
		     abs( start_nodes_ref - (start_node + pad_size) ) <= max_shift );
%  goodpoints = find( abs( start_nodes_ref - start_node ) <= max_shift );
  
  start_nodes_ref_all{ n } = start_nodes_ref( goodpoints );

  conv_matrix_all{ n } = conv_matrix( goodpoints, :);

  [xgrid,ygrid] = meshgrid( block_sizes, start_nodes_ref( goodpoints  ));
  final_nodes_ref = xgrid + ygrid - 1;
  % Disallow ending points that are outside the window (into the padding).
  badpoints = find( final_nodes_ref > (final_node_ref+pad_size) );

  % little check
  %s = find( scales == 1.0 );
  %[start_node final_node start_node_ref final_node_ref ]
  %[ min( start_nodes_ref(goodpoints))  min( final_nodes_ref(s,:)) final_node_ref+pad_size]
  %length( find( final_nodes_ref( s, : ) <= (final_node_ref + pad_size) ) )
  
  final_nodes_ref( badpoints ) = 2 * num_pixels_pad; % push to high number
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
DUMMY_NUMBER = -9999999999;
DP_matrix  = DUMMY_NUMBER * ones( num_pixels_pad, NUM_WINDOWS);
choice = zeros( num_pixels_pad, NUM_WINDOWS );

DP_matrix( :, NUM_WINDOWS+1) = 0.0;
tic
for n = NUM_WINDOWS:-1:1 
 
  fprintf('Dynamic programming for window: %d of %d\n',n, NUM_WINDOWS);

  start_nodes = start_nodes_ref_all{ n };
  final_nodes = final_nodes_ref_all{ n };
  conv_matrix = conv_matrix_all{n};

  for k = 1:length( start_nodes )

    in_range_points = find( final_nodes(k,:) < num_pixels_pad );
    in_DP_points = find(  DP_matrix( final_nodes(k,in_range_points)+1, n+1 ) ~= DUMMY_NUMBER  )';
    goodpoints = in_range_points( in_DP_points );
    if length( goodpoints ) > 0
      [y, i ]  = max( conv_matrix( k, goodpoints ) + ...
		      DP_matrix( final_nodes(k, goodpoints)+1, n+1)' );      
      DP_matrix( start_nodes(k) , n) = y;      
      choice( start_nodes(k), n) = final_nodes( k, goodpoints(i) );
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

