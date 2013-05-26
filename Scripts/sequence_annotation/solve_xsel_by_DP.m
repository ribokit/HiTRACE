function [ xsel_fit, best_score ] = solve_xsel_by_DP( I_data, alpha_ideal, sequence_at_bands, ideal_spacing, input_bounds, data_types )
% [ xsel_fit, best_score ] = solve_xsel_by_DP( I_data, alpha_ideal, sequence_at_bands, ideal_spacing, input_bounds, data_types );
%
% Check autoassign_notes_on_leastsquares_DP_rhiju.pdf  (checked into HiTRACE repository in same directory as this script)
% for description of the 'scorefunction'.
%
% (C) R. Das 2010-2013

if nargin == 0;  help( mfilename ); return; end;

START_POS = 0;
END_POS = 0;

% basic initialization
N = size( alpha_ideal, 1 );
[num_pixels, num_lanes] = size( I_data );
x = (1:num_pixels)';

if length( input_bounds ) >= 1; START_POS = round(input_bounds(1)); end;
if length( input_bounds ) >= 2; END_POS   = round(input_bounds(end)); end;

MIN_SEP = round(ideal_spacing/2); % this doesn't make sense -- optimal_SEP = ideal_spacing/3?
MAX_SEP = round(ideal_spacing*1.5);

SPACING_WEIGHT = 0.005; %0.1; % 0.1;
XSEL_GIVEN = 0;
if length( input_bounds ) == N; 
  XSEL_GIVEN = 1;
  spacings =  abs(input_bounds(1:end-1) - input_bounds(2:end));
  ideal_spacing = round( mean(spacings ) ); 
  MIN_SEP = min( spacings );
  MAX_SEP = max( spacings );
end;

% Define a model for peak separations
peak_width = ideal_spacing/6.0; % gaussian width.
beta_SEP = SPACING_WEIGHT * ones(N, 1); % weights for each spacing bonus
optimal_SEP = ideal_spacing * ones(N, 1); % ideal spacing for each spacing bonus
% band compression at G's.
for n = 1:(N-1);
  if ( sequence_at_bands( n ) == 'G' )
    beta_SEP( n+1 )    = beta_SEP( n+1 ) / 2;
    optimal_SEP( n+1 ) = optimal_SEP( n+1 ) / 3;
  end
end

% We need to normalize I_data, right?
%I_data = I_data/mean(mean(I_data));
%I_data = (max(I_data,0)/5).^0.5;
%I_data = max(I_data,0)/5;
window_size = ideal_spacing * 2;
I_data = window_normalize( I_data, window_size )/2;

% Set up basis functions.
% A gaussian. 
xmid = round(median(x));
ideal_gaussian = get_gaussian( x, xmid, peak_width );
g0 = sum(ideal_gaussian.^2);
g = zeros(MAX_SEP,1);
% Gaussian-gaussian correlation (needed for Ipred^2 calc below).
for i = 1:MAX_SEP
  g(i) = sum( ideal_gaussian .* get_gaussian(x, xmid+i, peak_width ) );
end

% Set up a 'peak bonus'
PEAK_WEIGHT = 1; % 1;
PEAK_MATRIX_WEIGHT = 0; % 1;
PEAK_SPREAD = round(peak_width/2); % number of pixels
ok_points = (PEAK_SPREAD+1):(num_pixels-PEAK_SPREAD);
peak_shifts = -PEAK_SPREAD:PEAK_SPREAD;
peak_scores = ones(1,num_lanes);
if length( data_types ) > 0
  for i = 1:num_lanes
    if strcmp(data_types{i},'nomod')
      %peak_scores(i) = 0;
    elseif ( strcmp(data_types{i},'ddTTP') || strcmp(data_types{i},'ddGTP') ...
	     || strcmp(data_types{i},'ddATP') || strcmp(data_types{i},'ddCTP') || strcmp(data_types{i},'ddUTP') )
      %peak_scores(i) = (num_lanes - 2) / 1.5;
    end
  end
end

tic
fprintf( 'Computing peak bonus matrix...\n');
peak_scores = peak_scores ./ length(peak_shifts);
peak_bonus_matrix = zeros( num_pixels, num_lanes );
for n = 1:num_lanes
  % this peak finder takes some time...
  peaks = getpeaks( I_data(:,n), 'THRESHOLD', min([ mean(I_data(:,n)), median(I_data(:,n)), max(I_data(:,n))*0.05 ]) );
  peaks = intersect( peaks, ok_points );
  for k = peaks
    peak_bonus_matrix( peak_shifts+k, n ) = peak_bonus_matrix( peak_shifts+k, n ) - peak_scores(n);
  end
end
peak_bonus = sum(peak_bonus_matrix, 2);
toc

%replace with FFT-based convolution
fprintf( 'Computing FFT-based convolution...\n');
n_fft_row = 2*num_pixels - 1;
gaussian_fft =  fftn( ideal_gaussian, [n_fft_row num_lanes] );
I_data_fft = fftn( I_data, [n_fft_row num_lanes]);
f = real( ifft2( gaussian_fft .* I_data_fft ) );
f = f( xmid + (0:num_pixels-1), : );


%figure; hold on; plot(I_data(:,1),'r');plot(I_data(:,2),'g');plot(I_data(:,3),'b');plot(I_data(:,5),'k');plot(peak_bonus,'k');
%figure; hold on; plot(f(:,1),'r');plot(f(:,2),'g');plot(f(:,3),'b');plot(f(:,5),'k');plot(peak_bonus,'k');

I_data_2 = sum( sum( I_data.^2 ) );
D = nan * ones( num_pixels, N );
I_pred_data_DP = nan * ones( num_pixels, N );
I_pred2_DP = nan * ones( num_pixels, N );
prev_pos_best = nan * ones( size(D) ); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize matrix.
% I_data^2 + I_pred^2 - 2 * I_data*I_pred
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start_vals = I_data_2  +   sum(alpha_ideal(1,:).^2) * g0  -  2 * alpha_ideal(1,:) * f(x,:)' + PEAK_WEIGHT*peak_bonus';

% first peak position.
if ( START_POS > 0 ) % only one value filled.
  D(START_POS,1) = start_vals( START_POS );  
else
  D(:,1) = start_vals;
end

I_pred_data_DP( :,1) = -2 * alpha_ideal(1,:) * f(x,:)';
I_pred2_DP( :,1) = sum(alpha_ideal(1,:).^2) * g0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fill matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% used for calculating peak^2 [for predicted profiles]
sum_n = sum(alpha_ideal.^2, 2);
% used for calculating peak-peak overlap [for predicted profiles]
overlap_n_nminus1 = sum(alpha_ideal .* [zeros(1,num_lanes); alpha_ideal(1:(end-1),:)], 2); 
D_overlap1 = sum_n * g0;
D_overlap2 = 2 * overlap_n_nminus1;

min_pos_raw = 0;
max_pos_raw = num_pixels - sum(optimal_SEP) + optimal_SEP(1);

CHECK_NUMERIC = 0;

for n = 2:N
  
  if ~XSEL_GIVEN, fprintf( 'Checking band position: %d of %d\n', n, N ); end;
  
  min_pos_raw = min_pos_raw + optimal_SEP(n-1);
  max_pos_raw = max_pos_raw + optimal_SEP(n);
    
  min_pos = floor(min_pos_raw);
  reasonable_pixel_range = [ min_pos : floor(max_pos_raw) ];  
  
  if XSEL_GIVEN % entire xsel is already provided
    reasonable_pixel_range = [ input_bounds(n) ];
  end
    
  prev_pos_min = max(reasonable_pixel_range - MAX_SEP, 1)';
  prev_pos_max = max(reasonable_pixel_range - MIN_SEP, 1)';
    
  D_datapred_peak_peakpred = nan * ones( 1, num_pixels );
  D_datapred_peak_peakpred( reasonable_pixel_range ) = ...
      - 2 * f(reasonable_pixel_range,:) * alpha_ideal(n,:)' ... % score for Idata*Ipred
      + PEAK_WEIGHT * num_lanes * peak_bonus(reasonable_pixel_range) ... % bonus for being at a peak location.
      + PEAK_MATRIX_WEIGHT * num_lanes * peak_bonus_matrix(reasonable_pixel_range,:) * alpha_ideal(n,:)'; % score for peak*Ipred
    
  % following could probably be sped up as a meshgrid calculation, or at least parallelized.
  for i = 1:length(reasonable_pixel_range)
    
    new_peak_pos = reasonable_pixel_range( i );

    % look at possible 'previous positions' -- this is in global coordinates
    prev_pos = [prev_pos_min(i) : prev_pos_max(i)];

    % check where prev_pos is not nan!!!!!!!
    prev_pos = prev_pos(  ~isnan( find( D(prev_pos,n-1) ) ) );

    spacings = new_peak_pos - prev_pos'; % these are spacings.

    D_test = nan * ones(1,num_pixels);
    D_test( prev_pos ) = ...
	D( prev_pos, n-1 ) ...
	+ (D_overlap1(n) + D_overlap2(n) * g( spacings )) ...
	+ beta_SEP(n) * num_lanes * ( spacings - optimal_SEP(n) ).^2  ...
	+ D_datapred_peak_peakpred( new_peak_pos );

    [ D(new_peak_pos,n),  best_idx ] = min( D_test(prev_pos) );
    prev_pos_best(new_peak_pos,n) = prev_pos( best_idx ); 
    
    % for numerical consistency check
    if CHECK_NUMERIC
      I_pred_data_DP( new_peak_pos, n ) =  I_pred_data_DP( prev_pos(best_idx), n-1) - 2 * f( new_peak_pos,:)*alpha_ideal(n,:)';
      I_pred2_DP( new_peak_pos, n ) =  I_pred2_DP( prev_pos(best_idx), n-1) + (D_overlap1(n) + D_overlap2(n) * g( spacings(best_idx) ));
    end
    
  end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Backtrack
xsel_fit = [];
[best_score, xsel_fit(N) ] = min( D(:,N) );
fprintf( 'Found best score: %8.1f\n', best_score );

if ( END_POS > 0 ); xsel_fit(N) = END_POS; end;

for n = (N-1) : -1 : 1
  xsel_fit(n) = prev_pos_best( xsel_fit(n+1), n+1 );
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% visual feedback & cross checks
I_pred = get_Ipred(alpha_ideal, xsel_fit,x,peak_width); 

if false, show_DP_matrix( D,xsel_fit); end;
if false, show_fit_matrix( I_data, I_pred, xsel_fit, alpha_ideal, peak_bonus_matrix ); end;
if CHECK_NUMERIC, do_numerical_checks( I_data, I_pred, alpha_ideal, xsel_fit, f, g, g0, I_data_2, I_pred2_DP, I_pred_data_DP ); end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_xsel( xsel, alpha_ideal )
hold on
for i = 1:size( alpha_ideal, 2 );
  for m = 1:length(xsel )
    if alpha_ideal(m,i) > 0.5
      plot( i, xsel(m), 'ro' );
    end
  end
end
hold off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% some leftover plotting funcitons -- might be useful in the future
function show_DP_matrix( D, xsel_fit); 
N = size(D,2);
D(:,1) = [];
D = D - min(min(D)) + 1;
%D(D>92891) = Inf;
D = (log(D) - 7.85) * 70;
%D = D - min(min(D));
%fprintf('%.3f\n',max(D(D<Inf)));
%D = ( D - 1000 ) / 80;

figure;hold on;image( D );colormap('Pink');colorbar();
for n = N : -1 : 1
  plot(n-1,xsel_fit(n),'ro','MarkerSize',7,'LineWidth',2);
end
set(gca,'YDir','reverse');
axis tight;
xlim([0.5 N-0.5]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function I_pred = get_Ipred( alpha_ideal, xsel, x, peak_width );

num_lanes = size( alpha_ideal, 2 );
num_peaks = length( xsel );
I_pred = zeros( length(x), num_lanes ) ;
for j = 1:num_lanes
  for i = 1:num_peaks;
    I_pred(:,j) = I_pred(:,j) + alpha_ideal(i,j) * get_gaussian( x, xsel(i), peak_width );
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function show_fit_matrix( I_data, I_pred, xsel, alpha_ideal, peak_bonus_matrix );

figure(2)
num_lanes = size( alpha_ideal, 2 );

subplot( 1, 3, 1 );
image( I_data * 100 );
plot_xsel( xsel, alpha_ideal );
title( 'data' );

subplot( 1, 3, 2 );
image( peak_bonus_matrix * -500 );
plot_xsel( xsel, alpha_ideal );
title( 'peak bonus matrix' );

subplot( 1, 3, 3 );
image( I_pred* 100 );
%plot_xsel( xsel, alpha_ideal );
title( 'fit' );

colormap( 1 - gray(100) );

%plot_xsel( xsel, alpha_ideal );

figure(3)  
Nplot = min( 4, num_lanes );
for n = 1:Nplot
  subplot( Nplot,1,n);
  plot( I_data(:,n), '.' ); hold on
  plot( I_pred(:,n) ); 
  for m = 1:length( xsel )
    plot( xsel(m)*[1 1], [0 I_pred(xsel(m),n)],'k' );
  end
  hold off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function do_numerical_checks( I_data, I_pred, alpha_ideal, xsel_fit, f,  g, g0, I_data_2, I_pred2_DP, I_pred_data_DP );

num_lanes = size( alpha_ideal, 2 );
N = length( xsel_fit );

alpha_f_sum = 0.0;
for j = 1:num_lanes
  alpha_f_sum = alpha_f_sum + sum( alpha_ideal(:,j) .* f( xsel_fit, j ) );
end
fprintf( 'Ipred*Idata ==> %8.1f (actual) or %8.1f (gauss) or %8.1f\n', sum( sum( I_data.*I_pred ) ), alpha_f_sum, -0.5*I_pred_data_DP( xsel_fit(N), N) );
  
% Check norm of Ipred... tests pairwise gaussian sum.
spacings = xsel_fit(2:end) - xsel_fit(1:end-1);
I_pred2_gauss = g0 * sum( sum( alpha_ideal.^2 ) )  +  2*  sum( g( spacings ).* sum(alpha_ideal(1:end-1,:) .* alpha_ideal(2:end,:),2 ) );
fprintf( 'Ipred*Ipred ==> %8.1f (actual) or %8.1f (gauss) or %8.1f (DP)\n', ...
	 sum( sum( I_pred .* I_pred ) ),  ...
	 I_pred2_gauss, ...
	 I_pred2_DP( xsel_fit(N), N ) ...
	 );

dev_DP = I_data_2 + I_pred2_gauss - 2 * alpha_f_sum;

fprintf( '(Ipred-Idata)^2 ==> %8.1f (actual) or %8.1f (gauss)\n', sum(sum( (I_data-I_pred ).^2 )), dev_DP );
