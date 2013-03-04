function [d_out,x_transform_all, anchor_nodes] = align_by_DP( d_all, align_blocks_in, penalizeStretchFactor, slack, maxShift, windowSize,  PLOT_STUFF );
% ALIGN_BY_DP: refine alignment by piece-wise-linear transform, optimizing correlation by dynamic programming
%
%  [d_out, x_transform_all, anchor_nodes] = align_by_DP( d_all, align_blocks_in, penalizeStretchFactor, slack, maxShift, windowSize,  PLOT_STUFF );
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
%  x_transform_all   = [advanced] matrix describing the local realignments
%  anchor_nodes = [advanced] window boundaries. 
%
% (C) R. Das, 2009-2011, 2013
%

d_out = [];
if nargin == 0;  help( mfilename ); return; end;

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
if ~exist( 'slack' ); slack = 50; end;
if ~exist( 'maxShift' ); maxShift = 200; end;
if ~exist( 'windowSize' ); windowSize = 500; end;
if ~exist( 'PLOT_STUFF' ); PLOT_STUFF = 1; end;

% need to double-check that the align_blocks are acceptable...
nlanes = size( d_all, 2 );
for j = 1:length( align_blocks_in )
  align_blocks_in{j} = min(max(align_blocks_in{j},1),nlanes );
end

d_out = d_all;

for j = 1:length( align_blocks_in )
  [d_align, x_transform_all, DP, choice, anchor_nodes ] = align_by_DP_block( d_out(:, align_blocks_in{j}), 1, penalizeStretchFactor, slack, maxShift, windowSize, PLOT_STUFF );
  d_out(:, align_blocks_in{j} )  = d_align;
end

if PLOT_STUFF
  colormap( 1 - gray(100 ) );
  subplot(1,2,1);
  scalefactor = 40 / mean(mean(d_all));
  image( scalefactor * d_all )
  show_anchor_nodes( anchor_nodes );
  %make_lines( [0:size(d_out,2)], 'k', 0.25 );
  
  subplot(1,2,2);
  image( scalefactor * d_out )
  %make_lines( [0:size(d_out,2)], 'k', 0.25 );
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [d_out,x_transform_all,DP,choice, anchor_nodes] = align_by_DP_block( d_all, refcol, penalizeStretchFactor, slack, maxShift, windowSize, PLOT_STUFF );

DP = [];
choice = [];

%d_ref = d_ref - mean( d_ref );
%d = d - mean( d );
which_lanes =  1:size( d_all, 2 );

num_pixels = size( d_all, 1 );

x_transform_all = zeros(  num_pixels, size( d_all, 2 ) );

% keep track of alignment 'anchor' nodes.
nodes = get_windows( windowSize, d_all(:,1) );
anchor_nodes = zeros( length( nodes )+1, size( d_all, 2) );

% parallelization! yea!
if parallelization_exists()
  %for i = which_lanes;
  parfor i = which_lanes;
    [x_transform_all(:,i), anchor_nodes(:,i)]    = align_by_DP_inner_loop( i, refcol, d_all, penalizeStretchFactor, slack, maxShift, windowSize, num_pixels );
  end
else
  for i = which_lanes;
    [x_transform_all(:,i), anchor_nodes(:,i)]    = align_by_DP_inner_loop( i, refcol, d_all, penalizeStretchFactor, slack, maxShift, windowSize, num_pixels );
  end
end
  
d_out = apply_transform( d_all, x_transform_all, PLOT_STUFF );
   
return;  
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
function  [x_transform, anchor_nodes]  = align_by_DP_inner_loop( i, refcol, d_all, penalizeStretchFactor, slack, maxShift, windowSize, num_pixels );
  
if ( i == refcol )
  x_transform = [1:num_pixels];
  [start_nodes, final_nodes] = get_windows( windowSize, d_all(:,i) );
  anchor_nodes = [start_nodes, final_nodes(end) ];
else
  fprintf( '-- Aligning %d of %d --\n',i,size(d_all,2));
  %tic
  
  d_ref = process_profile( d_all(:,refcol) );
  d_ali = process_profile( d_all(:,     i) );
    
  %plot( [d_ref d_ali] ); pause;
  
  [d_transform, x_transform, DP, choice, anchor_nodes] = refine_by_transforming( d_ref, d_ali, penalizeStretchFactor, slack, maxShift, windowSize );    
  %toc
end

x_transform = x_transform';



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ d_transform, x_transform, DP_matrix, choice, anchor_nodes ] = refine_by_transforming( d_ref_input, d_ali_input, penalizeStretchFactor, slack, max_shift, window_size );

% window_size = 100; % in pixels
% slack = 10; % Stretch/contract each window by this many pixels
% max_shift = 100; % how far can this window drift? 

[start_nodes_ref, final_nodes_ref, NUM_WINDOWS, num_pixels] = get_windows( window_size, d_ref_input );

DP_matrix = zeros( num_pixels, NUM_WINDOWS );

% Pad alignment data with zeros on either end.
pad_size = max_shift;
d_ali = zeros( pad_size + num_pixels + pad_size, 1 );
relevant_middle_pixels = pad_size + [1:num_pixels];
d_ali( relevant_middle_pixels ) = d_ali_input;

num_pixels_pad  = length( d_ali );

norm2_ali = sum( d_ali.*d_ali );

parts_that_count = d_ali * 0;
parts_that_count( relevant_middle_pixels ) = 1;

% Precompute all necessary correlation values.
%tic
start_nodes_ali_PAD_all = {};
final_nodes_ali_PAD_all = {};
diff2_matrix_all = {};

for n = 1: NUM_WINDOWS

  %fprintf( 'Correlation calculation for window: %d of %d\n',n, NUM_WINDOWS);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Window of profile from reference signal
  start_node_ref = start_nodes_ref( n );
  final_node_ref = final_nodes_ref( n );
  
  window_ref = [start_node_ref : final_node_ref];
  window_length = final_node_ref - start_node_ref + 1;

  min_window_length = max( window_length - slack, 2);
  max_window_length = window_length + slack;

  d_ref_sub  = d_ref_input( window_ref );

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Stretch this segment.
  % reference signal with different stretching scale factors
  block_sizes = [min_window_length:1:max_window_length];
  [d_ref_scaled, parts_that_count, scales] = get_d_ref_scaled( block_sizes, window_length, d_ref_sub );
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Carve out a window of the profile-alignment that hopefully
  % encapsulates a segment corresponding to the window of the reference
  % add pad_size because of zero padding at beginning of alignment reference.
  start_node_ali = pad_size + (start_node_ref - max_shift);
  final_node_ali = pad_size + (final_node_ref + max_shift);
  num_pixels_ali_window = final_node_ali - start_node_ali + 1;
  d_ali_window = d_ali( [start_node_ali:final_node_ali] );
    
  % Do convolution calculation. Need zero padding! 
  n_fft_row = max( length( d_ali_window ), max_window_length );
  n_fft_col = length( scales );

  % do the padding explicitly (fft should do it, but I think it might be symmetric?)
  d_ref_scaled     = pad_matrix( d_ref_scaled, n_fft_row );
  d_ali_window     = pad_matrix( d_ali_window, n_fft_row );
  parts_that_count = pad_matrix( parts_that_count, n_fft_row );
   
  % Following approximates:  [profile1 - profile2 ] ^2 
  diff2_matrix = get_diff2_matrix( d_ali_window, d_ref_scaled, parts_that_count, n_fft_row, n_fft_col );

  % should also scale intensities -- otherwise stretched windows will give a bonus -- they hit more pixels.
  diff2_matrix = diff2_matrix * diag( 1./scales ) ;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %  Start nodes and final nodes -- note that this retains max_shift offset!
  start_nodes_ali = start_node_ali - 1 + [1:length(d_ali_window)];
  
  % Disallow shifts greater than max_shift.
  goodpoints_start = find( abs( start_nodes_ali - (start_node_ref + max_shift) ) <= max_shift );
  start_nodes_ali = start_nodes_ali( goodpoints_start );
  
  start_nodes_ali_all{ n } = start_nodes_ali;
  diff2_matrix_all{ n } = diff2_matrix( goodpoints_start, :);

  [xgrid,ygrid] = meshgrid( block_sizes, start_nodes_ali );
  final_nodes_ali = xgrid + ygrid - 1;

  % Disallow ending points that are outside the window. 
  badpoints = find( abs( final_nodes_ali  - (final_node_ref + max_shift) ) > max_shift );
  final_nodes_ali( badpoints ) = NaN;
  final_nodes_ali_all{ n } = final_nodes_ali;

  % plots to check crazy fft convolutions above
  if ( n == 0 ) %>= NUM_WINDOWS-2 )
    check_scale = 1.0;
    k = find( scales == check_scale);
    d_ref_scaled_show = d_ref_scaled(:,k );
    subplot(5,1,1); plot( d_ref_scaled_show ) ;  title( 'd_ref_scaled [REF] -- scale 1.0')
    xlim([0 n_fft_row])    
    subplot(5,1,2); plot( d_ali_window ); title( 'd-ali-window')
    xlim([0 n_fft_row])
    subplot(5,1,3); plot( diff2_matrix(:,k) );
    xlim([0 n_fft_row])
    subplot(5,1,4); plot( norm2_over_window(:,find( scales == 1.0) ) );
    xlim([0 n_fft_row])
    subplot(5,1,3); 
    for m = 1: length( d_ali_window )
      diff2_test( m ) = sum( (circshift(d_ref_scaled_show,m) - circshift( parts_that_count(:,k),m ).*d_ali_window).^2 );
    end
    
    hold on
    plot( diff2_test * check_scale,'r' ); 
    plot( goodpoints_start, diff2_test( goodpoints_start ) * (1/check_scale),'g' ); 
    hold off

    subplot(5,1,5);
    plot( min( diff2_matrix' ) );
    [dummy, idx_j] = min( diff2_matrix( goodpoints_start,:) );
    [minval,idx ]=  min( dummy )
    
    %start_nodes_ali( idx_j( idx ) )
    %image( 0.05 * diff2_matrix( 1:1000,:)'  );
    %plot( d_ref_scaled )

    pause;
  end

end
%toc

% with the conventions I chose, easier to work backward, last
% window to first.
DP_matrix  = NaN * ones( num_pixels_pad+1, NUM_WINDOWS);
choice          = zeros( num_pixels_pad, NUM_WINDOWS );
best_block_stretch = zeros( num_pixels_pad, NUM_WINDOWS );

% stretch penalty. Need to make it comparable to scale of profile convolution score
PENALIZE_STRETCH = penalizeStretchFactor * ( std( d_ref_input ) )^2 * ( slack/window_size)^2;

DP_matrix( :, NUM_WINDOWS+1) = 0.0;
%tic
for n = NUM_WINDOWS:-1:1 


  %fprintf('Dynamic programming for window: %d of %d\n',n, NUM_WINDOWS);

  start_nodes = start_nodes_ali_all{ n };
  final_nodes = final_nodes_ali_all{ n };
  diff2_matrix = diff2_matrix_all{n};

  for k = 1:length( start_nodes )

    in_range_points = find( ~isnan( final_nodes(k,:) ) );
    in_DP_points = find(  ~isnan( DP_matrix( final_nodes(k,in_range_points)+1, n+1 ) )  )';
    goodpoints = in_range_points( in_DP_points );
    
    if length( goodpoints ) > 0
      
      prev_score = DP_matrix( final_nodes(k, goodpoints)+1, n+1)';
      convolution_score = diff2_matrix( k, goodpoints );
      % New -- prevent transforming if there's not much information.
      stretch_score = 0 * prev_score;
      potential_block_stretches = ( final_nodes(k, goodpoints) - start_nodes(k) ) - ( final_nodes_ref(n) - start_nodes_ref(n) );
      if ( n < NUM_WINDOWS )
	%stretch_score = PENALIZE_STRETCH * abs(potential_block_stretches - best_block_stretch( final_nodes(k, goodpoints)+1, n+1 )');
	stretch_score = PENALIZE_STRETCH * (potential_block_stretches - best_block_stretch( final_nodes(k, goodpoints)+1, n+1 )').^2;
      end
      
      [y, i ]  = min( prev_score + convolution_score + stretch_score );      
      DP_matrix( start_nodes(k) , n) = y;      
      choice( start_nodes(k), n) = final_nodes( k, goodpoints(i) );
      best_block_stretch( start_nodes(k), n ) = potential_block_stretches( i );
    end


  end

  if( n == 0 )
    clf; 
    plot( DP_matrix(:,n) ); 
    gp = find( ~isnan( DP_matrix(:,n) ) );
    [dummy,idx] = min( DP_matrix( gp, n ) );
    gp( idx )
    %pause; 
  end;

end
%toc

%clf
%plot( DP_matrix )
%pause;

% Backtrack.
[y, node] = min( DP_matrix(:,1) );

start_nodes = [];
final_nodes = [];
for n = 1:NUM_WINDOWS
  start_nodes(n) = node;  
  stretches(n) =   best_block_stretch( node, n );
  node = choice( node, n )+1;
  final_nodes(n) = node-1;
end


% Remove max_shift offset!
start_nodes = start_nodes - max_shift + 1;
final_nodes = final_nodes - max_shift + 1;

%stretches
%start_nodes
%final_nodes

nodes = [ start_nodes final_nodes(end) ];
nodes_ref = [ start_nodes_ref final_nodes_ref(end)];
%[ nodes; nodes_ref; nodes-nodes_ref]

x_transform = interp1( nodes, nodes_ref, [1:num_pixels], 'linear','extrap');

if ( length( unique( x_transform )  ) ~= length( x_transform  ) ) % something weird.
  x_transform;
  x_transform = 1:num_pixels;
end
d_transform = interp1( x_transform, d_ali_input( 1: num_pixels), 1:num_pixels, 'linear',0.0);

anchor_nodes = nodes'; %useful output

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [start_nodes_ref, final_nodes_ref, NUM_WINDOWS, num_pixels] = get_windows( window_size, d_ref_input );

num_pixels = length( d_ref_input );

NUM_WINDOWS = floor( num_pixels / window_size ) + 1;

% Make sure last window is at least 1 pixel big.
if ( num_pixels - (NUM_WINDOWS-1)*window_size < 2 )
  NUM_WINDOWS = NUM_WINDOWS - 1;
end

for n = 1: NUM_WINDOWS
  start_nodes_ref(n) = (n-1)*window_size + 1;
  final_nodes_ref(n) = min( n*window_size, num_pixels );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  d =  pad_matrix( d, n_fft_row );
d = [d; zeros( n_fft_row - size(d,1), size(d,2) ) ];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = process_profile( d );

d = quick_norm( smooth(d) );

% put on scale of 0 to 100.
d = 20 * (max( d, 0) );

% remove baseline.
d = d - mode( round(smooth(d)) );

% get interquartile range, and use it to cap outliers.  
d = cap_outliers( d );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  d =  cap_outliers( d )
d_sort = sort( d );
N = length(d);
q3 = d_sort( round( N * 0.75 ) );
q1 = d_sort( round( N * 0.25) );
cutoff = 5*(q3-q1)  + q3;
d( find( d > cutoff ) ) = cutoff;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function   [d_ref_scaled, parts_that_count, scales] = get_d_ref_scaled( block_sizes, window_length, d_ref_sub );
max_window_length = max( block_sizes );

d_ref_scaled = zeros( max_window_length, length( block_sizes ) );
parts_that_count = zeros( max_window_length, length( block_sizes ) );
interp_vectors = zeros( max_window_length,  length( block_sizes ) );

for m = 1: length( block_sizes )
  incr = (window_length-1)/(block_sizes(m)-1);
  interp_vectors( [1:block_sizes(m)], m) = [1:incr:window_length]';
  parts_that_count( [1:block_sizes(m)], m ) = 1.0;
  scales(m) = 1 / incr;
end
d_ref_scaled = interp1( [1:window_length], d_ref_sub, interp_vectors, 'linear', 0.0 );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function   diff2_matrix = get_diff2_matrix( d_ali_window, d_ref_scaled, parts_that_count, n_fft_row, n_fft_col );

dot_product_matrix = real(ifft2(fftn( d_ali_window, [n_fft_row n_fft_col]) ...
				.* fftn(d_ref_scaled(end:-1:1,:), [n_fft_row n_fft_col])));
norm2_over_window  = real(ifft2(fftn( d_ali_window.^2, [n_fft_row n_fft_col]) ...
				.* fftn( parts_that_count(end:-1:1,:), [n_fft_row n_fft_col])));
norm2_ref = repmat( sum( d_ref_scaled.^2 ), n_fft_row, 1);

diff2_matrix = norm2_over_window - 2 * dot_product_matrix + norm2_ref;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function show_anchor_nodes( anchor_nodes );

for i = 1:size( anchor_nodes, 1 )
  hold on
  plot( [1:size(anchor_nodes,2)], anchor_nodes(i,:), 'r-' );
end
hold off;