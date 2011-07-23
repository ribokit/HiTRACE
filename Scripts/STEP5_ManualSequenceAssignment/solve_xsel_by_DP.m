function [ xsel_fit, D ] = solve_xsel_by_DP( I_data, alpha_ideal, sequence_at_bands, ideal_spacing, input_bounds, PLOT_STUFF );
% [ xsel_fit, D ] = solve_xsel_by_DP( I_data, alpha_ideal, sequence_at_bands, ideal_spacing, input_bounds, PLOT_STUFF );

if ~exist( 'PLOT_STUFF' ); PLOT_STUFF = 1; end;
if ~exist( 'ideal_spacing'); ideal_spacing = 24; end;
if ~exist( 'input_bounds' ) input_bounds = []; end;

START_POS = 0;
if length( input_bounds ) >= 1; START_POS = round(input_bounds(1)); end;

END_POS = 0;
if length( input_bounds ) >= 2; END_POS = round(input_bounds(end)); end;

% basic initialization
width = ideal_spacing/4.0; % gaussian width.
N = size( alpha_ideal, 1 );
num_pixels = size( I_data, 1 );
num_lanes = size( I_data, 2 );
x = [1:num_pixels]';

% Define a model for peak separations
MIN_SEP = round(ideal_spacing/3);
MAX_SEP = ideal_spacing*2;
beta_SEP = 0.1*ones( N, 1 );
optimal_SEP = ideal_spacing * ones(N, 1);

% band compression at G's.
for n = 1:(N-1);
  if ( sequence_at_bands( n ) == 'G' )
    beta_SEP( n+1 ) = beta_SEP( n+1 ) /2;
    optimal_SEP( n+1 ) = optimal_SEP( n+1 ) / 3.0;
  end
end

% A gaussian. 
xmid = round(median(x));
ideal_gaussian = get_gaussian( x, xmid, width );
g0 = sum(ideal_gaussian.^2);

% Gaussian-gaussian correlation (needed for Ipred^2 calc below).
for i = 1:MAX_SEP
  g(i) = sum( ideal_gaussian .* get_gaussian(x, xmid+i, width ) );
end

% Gaussian-data correlation
%for i = x'
%  gaussian = repmat( get_gaussian(x, i, width ), [1 num_lanes] );
%  f( i, : ) = sum( I_data .* gaussian );
%end

%replace with FFT-based convolution
n_fft_row = 2*num_pixels - 1;
n_fft_col = num_lanes;
gaussian_fft =  fftn( ideal_gaussian, [n_fft_row n_fft_col ] );
I_data_fft = fftn( I_data, [n_fft_row n_fft_col]);
f = real( ifft2( gaussian_fft .* I_data_fft ) );
f = f( xmid + [0:num_pixels-1], : );

% Checks ... see sample code at end of script -- commented out

% Set up a 'peak bonus'
peak_bonus = zeros( 1, num_pixels );
PEAK_SPREAD = 1;
for n = 1:num_lanes
  peaks = localMaximum( I_data(:,n), ideal_spacing );
  ok_points = [(PEAK_SPREAD+1):(num_pixels-PEAK_SPREAD)];
  peaks = intersect( peaks,  ok_points );
  peak_shifts = [-PEAK_SPREAD:PEAK_SPREAD];
  for k = peak_shifts
    peak_bonus( peaks+k ) = peak_bonus( peaks+k ) - 1./length(peak_shifts);
  end
end
PEAK_WEIGHT = 1.0;

%BIG_NUMBER = 1e10;
I_data_2 = sum( sum( I_data.^2 ) );
BIG_NUMBER = 4*I_data_2 + sum( beta_SEP ) * MAX_SEP^2;
D = BIG_NUMBER * ones( num_pixels, N );

prev_pos_best = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize matrix.
% I_data^2 + I_pred^2 - 2 * I_data*I_pred
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start_vals = I_data_2  +   sum(alpha_ideal(1,:).^2) * g0  -  2 * alpha_ideal(1,:) * f(x,:)' + PEAK_WEIGHT*peak_bonus;

if ( START_POS > 0 )
  START_POS
  D(START_POS,1) = start_vals( START_POS );  
else
  D(:,1) = start_vals;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fill matrix
for n = 2:N
  fprintf( 'Filling %d of %d\n',n,N);

  min_pos = round( 0 + sum( 0.8*optimal_SEP( 1:(n-1) ) ) );
  max_pos = round( num_pixels - sum( 0.8*optimal_SEP( (n+1):end ) ) );

  reasonable_pixel_range = [min_pos:max_pos]; 
  %reasonable_pixel_range = [n:num_pixels];
  
  for i = reasonable_pixel_range

    D_test = 10*BIG_NUMBER*ones(1,num_pixels);

    % look at possible 'previous positions'
    prev_pos_min = max(i-MAX_SEP,1);
    prev_pos_max = max(i-MIN_SEP,1);
    prev_pos = [prev_pos_min:prev_pos_max];
    
    blah = 1;

    sum_n = sum(alpha_ideal(n,:).^2);
    overlap_n_nminus1 = sum(alpha_ideal(n,:) .* alpha_ideal(n-1,:));

    D_overlap = ...
	+ sum_n * g0 + ... % score for Ipred^2
	+ ( 2 * overlap_n_nminus1 * g(i - prev_pos)' ) ... % score for Ipred^2
	- 2 * sum( alpha_ideal(n, :) .* f( i, : ) ); % score for Idata*Ipred       

    D_test( prev_pos )  = ...	
	D( prev_pos, n-1) ...
	+ D_overlap/num_lanes ...
	+ beta_SEP(n) * ( i - prev_pos' - optimal_SEP(n) ).^2 ... % Score to maintain 'optimal' separation
	+ PEAK_WEIGHT * peak_bonus( i ); % bonus for being at a peak location.
    
    [ D(i,n),  prev_pos_best(i,n) ] = min( D_test );
    
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Backtrack
xsel_fit = [];
image( prev_pos_best )
image( D/1000)
[Dmin, xsel_fit(N) ] = min( D(:,N) );
if ( END_POS > 0 ) xsel_fit(N) = END_POS; end;

for n = (N-1) : -1 : 1
  xsel_fit(n) = prev_pos_best( xsel_fit(n+1), n+1 );  
end










% Checks ... 
%if exist( 'xsel' )
%
%  for j = 1:num_lanes
%    for i = 1:length( xsel );
%      I_pred(:,j) = I_pred(:,j) + alpha_ideal(i,j) * get_gaussian( x, xsel(i), width );
%    end
%  end
%  
%  % Check dot product
%  subplot(2,1,2)
%  plot( f );
%  
%  alpha_f_sum = 0.0;
%  for j = 1:num_lanes
%    alpha_f_sum = alpha_f_sum + sum( alpha_ideal(:,j) .* f( xsel, j ) );
%  end
%
%  fprintf( 'Ipred*Idata ==> %8.1f (actual) or %8.1f (gauss)\n', sum( sum( I_data.*I_pred ) ), alpha_f_sum );
%  
%  % Check norm of Ipred... tests pairwise gaussian sum.
%  delta = xsel(2:end) - xsel(1:end-1);
%% fprintf( 'Ipred*Ipred ==> %8.1f (actual) or %8.1f (gauss)\n', sum( sum( I_pred .* I_pred ) ),  g0 * sum( sum( alpha_ideal.^2 ) )  +  2*  sum( g( delta ).* sum(alpha_ideal(1:end-1,:) .* alpha_ideal(2:end,:),2 )' ) );
%
%end

