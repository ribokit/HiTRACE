function [data_align, x_realign, d1] = align_to_first_ver3( data, PLOT_STUFF, refcol )
%
% ALIGN_TO_FIRST_VER3:  (linear-time) alignment of a matrix of electropheretic traces to first trace
%
%   [data_align, x_realign, d1] = align_to_first_ver3( data, PLOT_STUFF, refcol )
%
%  Optimizes correlation coefficient by grid search + Fast Fourier Transform
%
% (C) R. Das & S.R. Yoon, 2009-2011
%

if nargin == 0;  help( mfilename ); return; end;

FULL_SIGNAL_WINDOW = 1000;

if ~exist('refcol', 'var');  refcol = 1; end

num_capillaries = size( data,2 );

if ~exist('PLOT_STUFF','var');  PLOT_STUFF = 0; end

%d_ref = baseline_subtract( data(:,refcol) );
d_ref = data(:,refcol);
[minbin_ref, middlebin_ref, maxbin_ref] = get_signal_bins( d_ref, FULL_SIGNAL_WINDOW );
%d1 = extract_profile( dezinger(d_ref) , minbin_ref, maxbin_ref );
d1 = extract_profile( d_ref , minbin_ref, maxbin_ref );
%d1 = d1 - smooth( d1, 100);

numpts_d_ref = length( d_ref );
x = [ 1: numpts_d_ref ];
x_realign = zeros(numpts_d_ref,num_capillaries);
data_align = zeros(numpts_d_ref,num_capillaries);
numpts = length(d1);
max_shift = 1000; min_shift = -1 * max_shift;
%scales = [0.93:0.005:1.08];
scales = [0.8:0.0025:1.2];
%scales = [1.0:0.05:1.4];

fprintf(' \n'); revStr = ' ';

for n = 1: num_capillaries
    %d     = baseline_subtract( data(:,n) );
    d     = data(:,n);
    
    [minbin, middlebin, maxbin ] = get_signal_bins( d, FULL_SIGNAL_WINDOW );
    %d2 = extract_profile( dezinger(d), minbin, maxbin );
    d2 = extract_profile( d, minbin, maxbin );
    %d2 = d2 - smooth( d2, 100);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    revStr = lprintf(revStr, ['Calculating correlation ', num2str(n), ' of ', num2str(num_capillaries), ' ... \n'], 2);
    
    % Do grid search
    [best_scale, best_shift] = find_best_scale_shift_COARSE( scales, min_shift, max_shift, d1, d2 );
    %[best_scale, best_shift] = find_best_scale_shift( scales, min_shift, max_shift, d1, d2 );
    
    new_x = best_scale * (x-minbin_ref+1 - best_shift)  + minbin - 1  ;
    da = interp1( x, d, new_x, 'linear',0);
    x_realign(:,n) = new_x;
    %da = baseline_subtract( da );
    data_align(:,n) = da;
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  [minbin, middlebin, maxbin ] =get_signal_bins_OLD( d, FULL_SIGNAL_WINDOW )
%Find continuous chunk with the most signal.

d = baseline_subtract( d );
d = d - smooth(d, 100 );

%numpts = length( d ); %sryoon

d_sum = cumsum( abs(d) );
d_sum2 = circshift( d_sum, FULL_SIGNAL_WINDOW );
tot_signal = d_sum - d_sum2;

[dummy, maxbin] = max( tot_signal((1+FULL_SIGNAL_WINDOW):end) );
maxbin = maxbin + FULL_SIGNAL_WINDOW;
minbin = maxbin - FULL_SIGNAL_WINDOW + 1;
middlebin = floor( 0.75*maxbin + 0.25*minbin  );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [minbin, middlebin, maxbin ] =get_signal_bins( d, FULL_SIGNAL_WINDOW )

minbin = 1;
maxbin = length( d );
middlebin = floor( 0.75*maxbin + 0.25*minbin  );

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d1 = extract_profile( d, minbin, maxbin )

%d1 = d( minbin:maxbin) - mean( d(maxbin+[1:100]) );
%d1 = peak_detect( d );
d1 = d;

d_sort = sort(d);
pctile_cutoff = 0.05;
d_min = d_sort( round( pctile_cutoff * length(d1)) );
d_max = d_sort( round( (1-pctile_cutoff) * length(d1)) );

d1 = max(min( d1, d_max),d_min);
%d1 = max( d1, 0 ).^0.25;

return;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tried following for speed increase -- no acceleration though!
function   [best_scale, best_shift] = find_best_scale_shift_COARSE( ...
    scales, min_shift, max_shift, d1, d2 )

bin_size = 20;

len_d1 = length(d1);

% probably a faster way to do this.
for i = 1: floor(len_d1 / bin_size)
    minbin = bin_size* ( i - 1 ) + 1;
    maxbin = min( bin_size* ( i ) + 1, len_d1 );
    d1_coarse( i ) = mean( d1( minbin:maxbin) );
    d2_coarse( i ) = mean( d2( minbin:maxbin) );
end

[ best_scale, best_shift ] = find_best_scale_shift( ...
    scales, ...
    min_shift/bin_size, max_shift/bin_size,...
    d1_coarse', d2_coarse');

best_shift = best_shift * bin_size;

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function   [best_scale, best_shift] = find_best_scale_shift( scales, ...
    min_shift, max_shift, d1, d2 )

numscales = length( scales );

len_d1 = length(d1);
n_fft_row = 2*len_d1 - 1;
n_fft_col = length(scales);
x = (1:len_d1)';
d2x=interp1(x, d2, x*scales, 'linear',d2(end));

% subtract out mean -- we want a correlation coefficient.
d1 = d1 - mean( d1 );
for k = 1:numscales
    d2x(:,k) = d2x(:,k) - mean( d2x(:,k));
    %norm_factor2( k ) = sqrt( sum( d2x(:,k).^2 ) );
end

% Necessary for correlation coefficient.
norm_factor1 = sqrt( sum( d1.^2 ) );
norm_factor2 = sqrt( sum( d2x.^2 ) );

% Sungroh Yoon's trick -- faster than matlab's conv2.
conv_matrix = real(ifft2(fftn(d1,[n_fft_row n_fft_col]) ...
    .* fftn(d2x(end:-1:1,:), [n_fft_row n_fft_col])));
%conv_matrix = ones(n_fft_row, n_fft_col);

shifts_all = [ -(len_d1-1):(len_d1-1) ];
goodpoints = find( shifts_all >= min_shift & shifts_all <= max_shift);
conv_matrix = conv_matrix(goodpoints,:);
shifts = shifts_all( goodpoints);

% It turns out that these normalization factors ( bringing
% corrcoeff in range 0 to 1) are really important.
corrcoeff_vs_scale = max( conv_matrix )./( norm_factor1 * norm_factor2);
% Put in a weak correction to keep scale factor close to 1.0.
%corrcoeff_vs_scale = corrcoeff_vs_scale - 10 * ( scales - 1.1).^2;
%plot( scales, corrcoeff_vs_scale );
%pause;

[ dummy, best_scale_index ] = max( corrcoeff_vs_scale );
[ dummy, best_shift_index ] = max( conv_matrix(:,best_scale_index) );
best_scale = scales( best_scale_index);
best_shift = shifts( best_shift_index );

