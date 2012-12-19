function [ xsel_fit, D, msg ] = solve_xsel_by_DP( I_data, alpha_ideal, sequence_at_bands, ideal_spacing, input_bounds )
% [ xsel_fit, D ] = solve_xsel_by_DP( I_data, alpha_ideal, sequence_at_bands, ideal_spacing, input_bounds );

START_POS = 0;
if length( input_bounds ) >= 1; START_POS = round(input_bounds(1)); end;

END_POS = 0;
if length( input_bounds ) >= 2; END_POS = round(input_bounds(end)); end;

% basic initialization
width = ideal_spacing/6.0; % gaussian width.
N = size( alpha_ideal, 1 );
[num_pixels, num_lanes] = size( I_data );
x = (1:num_pixels)';
msg = [];

% Define a model for peak separations
MIN_SEP = round(ideal_spacing/2);
MAX_SEP = round(ideal_spacing*1.5);
beta_SEP = 0.1 * ones(N, 1);
optimal_SEP = ideal_spacing * ones(N, 1);

% band compression at G's.
for n = 1:(N-1);
  if ( sequence_at_bands( n ) == 'G' )
    beta_SEP( n+1 ) = beta_SEP( n+1 ) / 2;
    optimal_SEP( n+1 ) = optimal_SEP( n+1 ) / 3;
  end
end

%tmp = sum(alpha_ideal, 2);
%optimal_SEP(tmp == 0) = optimal_SEP(tmp == 0) / 1.1;

% A gaussian. 
xmid = round(median(x));
ideal_gaussian = get_gaussian( x, xmid, width );
g0 = sum(ideal_gaussian.^2);
g = zeros(MAX_SEP,1);
% Gaussian-gaussian correlation (needed for Ipred^2 calc below).
for i = 1:MAX_SEP
  g(i) = sum( ideal_gaussian .* get_gaussian(x, xmid+i, width ) );
end

% Gaussian-data correlation
%for i = x'
%  gaussian = repmat( get_gaussian(x, i, width ), [1 num_lanes] );
%  f( i, : ) = sum( I_data .* gaussian );
%end

% Checks ... see sample code at end of script -- commented out

% Set up a 'peak bonus'
PEAK_SPREAD = 1;
PEAK_WEIGHT = 1.0;

ok_points = (PEAK_SPREAD+1):(num_pixels-PEAK_SPREAD);
peak_shifts = -PEAK_SPREAD:PEAK_SPREAD;
peak_scores = ones(1,num_lanes);
peak_scores(end-1) = 0;
peak_scores(end) = (num_lanes - 2) / 2;
peak_scores = peak_scores ./ length(peak_shifts);
peak_bonus_matrix = zeros( num_pixels, num_lanes );

for n = 1:num_lanes

    peaks = localMaximum( I_data(:,n), ideal_spacing/2 );
    peaks = peaks(I_data(peaks,n) > mean(I_data(peaks,n)/2));
    peaks = intersect( peaks,  ok_points );
    
    peaks_mat{n} = peaks';
    for k = peak_shifts
        peak_bonus_matrix( peaks+k, n ) = peak_bonus_matrix( peaks+k, n ) - peak_scores(n);
    end

end
peak_bonus = sum(peak_bonus_matrix, 2);

%replace with FFT-based convolution
n_fft_row = 2*num_pixels - 1;
gaussian_fft =  fftn( ideal_gaussian, [n_fft_row num_lanes] );
I_data_fft = fftn( I_data, [n_fft_row num_lanes]);
f = real( ifft2( gaussian_fft .* I_data_fft ) );
f = f( xmid + [0:num_pixels-1], : );

%figure; hold on; plot(I_data(:,1),'r');plot(I_data(:,2),'g');plot(I_data(:,3),'b');plot(I_data(:,5),'k');plot(peak_bonus,'k');
%figure; hold on; plot(f(:,1),'r');plot(f(:,2),'g');plot(f(:,3),'b');plot(f(:,5),'k');plot(peak_bonus,'k');

%BIG_NUMBER = 1e10;
I_data_2 = sum( sum( I_data.^2 ) );
BIG_NUMBER = 4*I_data_2 + sum( beta_SEP ) * MAX_SEP^2;
D = BIG_NUMBER * ones( num_pixels, N );
prev_pos_best = zeros(size(D));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize matrix.
% I_data^2 + I_pred^2 - 2 * I_data*I_pred
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start_vals = I_data_2  +   sum(alpha_ideal(1,:).^2) * g0  -  2 * alpha_ideal(1,:) * f(x,:)' + PEAK_WEIGHT*peak_bonus';

if ( START_POS > 0 )
  D(START_POS,1) = start_vals( START_POS );  
else
  D(:,1) = start_vals;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fill matrix
sum_n = sum(alpha_ideal.^2, 2);
overlap_n_nminus1 = sum(alpha_ideal .* [zeros(1,num_lanes); alpha_ideal(1:(end-1),:)], 2);

D_overlap1 = sum_n * g0;
D_overlap2 = 2 * overlap_n_nminus1;

D_test = 10 * BIG_NUMBER * ones(1,num_pixels);

min_pos_raw = 0;
max_pos_raw = num_pixels - sum(optimal_SEP) + optimal_SEP(1);
for n = 2:N
    
    min_pos_raw = min_pos_raw + optimal_SEP(n-1);
    max_pos_raw = max_pos_raw + optimal_SEP(n);
    
    min_pos = floor(min_pos_raw);
    reasonable_pixel_range = min_pos:floor(max_pos_raw);
    
    prev_pos_min = max(reasonable_pixel_range - MAX_SEP, 1)';
    prev_pos_max = max(reasonable_pixel_range - MIN_SEP, 1)';
    
    D_datapred_peak_peakpred = ...
        - 2 / num_lanes * f(reasonable_pixel_range,:) * alpha_ideal(n,:)' ... % score for Idata*Ipred
        + PEAK_WEIGHT * peak_bonus(reasonable_pixel_range) ... % bonus for being at a peak location.
        + num_lanes * peak_bonus_matrix(reasonable_pixel_range,:) * alpha_ideal(n,:)'; % score for peak*Ipred
    
    reasonable_pixel_range = reasonable_pixel_range - min_pos + 1;
    
    for i = reasonable_pixel_range

        rawindex = i + min_pos - 1;
        
        % look at possible 'previous positions'
        prev_pos = prev_pos_min(i):prev_pos_max(i);
        tmpindex = rawindex - prev_pos';
        
        %D_overlap = ...
        %    + D_overlap1(n) ... % score for Ipred^2
        %    + ( D_overlap2(n) * g(tmpindex) ); % score for Ipred^2
        
        %D_barrier = beta_SEP(n) * ( tmpindex - optimal_SEP(n) ).^2; % Score to maintain 'optimal' separation
        
%         if false
%             if n == 42 && D_datapred_peak_peakpred(i) < -1
%                 fprintf('--i:%d\n',i);
%                 fprintf('rawindex:%d\n',rawindex);
%                 fprintf('prev_pos:%d\n',prev_pos);
%                 fprintf('tmpindex:%d\n',tmpindex);
%                 fprintf('D_overlap:%f\n',D_overlap1(n) + D_overlap2(n) * g(tmpindex) / num_lanes);
%                 fprintf('sep:%f\n',beta_SEP(n) / num_lanes *( tmpindex - optimal_SEP(n) ).^2);
%                 fprintf('D_bonus2:%f\n',D_datapred_peak_peakpred(i));
%                 fprintf('D_overlap1:%f\n',D_overlap1(n));
%                 fprintf('D_overlap2:%f\n',D_overlap2(n));
%                 fprintf('D:%f\n',D( prev_pos, n-1 ));
%             end
%         end
                
        D_test(prev_pos) = ...
            D( prev_pos, n-1 ) ...
            + (D_overlap1(n) + D_overlap2(n) * g(tmpindex)) / num_lanes ...
            + beta_SEP(n) * ( tmpindex - optimal_SEP(n) ).^2 / num_lanes ...
            + D_datapred_peak_peakpred(i);
        
%         if false
%             
%             ttttt = num2cell(prev_pos);
%             for tttt = 2:length(prev_pos)-1
%                 ttttt{tttt} = '';
%             end
%             
%             if n == 39 && rawindex == 507
%                 figure;
%                 colormap([0.7 0.5 0.5; 0.5 0.5 0.7; 0.9 0.8 0.1; 0.6 0.7 0.7; 0.9 0.5 0.5; 0.55 0.55 0.55; 0.3 0.3 0.3]);
%                 subplot(1,4,1), barsh = barh([zeros(length(tmpindex),5), D(prev_pos,n-1), zeros(length(tmpindex),1)],'stacked'); set(gca,'YDir','reverse','XTickLabel',[],'YTickLabel',ttttt);
%                 axesh = get(barsh(1), 'Parent'); set(axesh, 'Xlim', [min(D(prev_pos,n-1))-2 max(D(prev_pos,n-1))+2]);
%                 subplot(1,4,2), barh([(D_overlap1(n) + D_overlap2(n) * g(tmpindex)) * 40, beta_SEP(n) * ( tmpindex - optimal_SEP(n) ).^2 / num_lanes,zeros(length(tmpindex),5)],'stacked'); set(gca,'YDir','reverse','XTickLabel',[],'YTickLabel',[]);
%                 subplot(1,4,3), barh(repmat([0,0,2 / num_lanes * f(rawindex,:) * alpha_ideal(n,:)', -PEAK_WEIGHT * peak_bonus(rawindex), -num_lanes * peak_bonus_matrix(rawindex,:) * alpha_ideal(n,:)',0,0] ...
%                     ,length(tmpindex),1),'stacked'); set(gca,'YDir','reverse','XTickLabel',[],'YTickLabel',[]);
%                 subplot(1,4,4), barsh = barh([zeros(length(tmpindex),6), D_test(prev_pos)'],'stacked'); set(gca,'YDir','reverse','XTickLabel',[],'YTickLabel',ttttt);
%                 axesh = get(barsh(1), 'Parent'); set(axesh, 'Xlim', [min(D_test(prev_pos))-2 max(D_test(prev_pos))+2]);
%             end
%             
%             if n == 49 && rawindex == 843
%                 figure;
%                 colormap([0.7 0.5 0.5; 0.5 0.5 0.7; 0.9 0.8 0.1; 0.6 0.7 0.7; 0.9 0.5 0.5; 0.55 0.55 0.55; 0.3 0.3 0.3]);
%                 subplot(1,4,1), barsh = barh([zeros(length(tmpindex),5), D(prev_pos,n-1), zeros(length(tmpindex),1)],'stacked'); set(gca,'YDir','reverse','XTickLabel',[],'YTickLabel',ttttt);
%                 axesh = get(barsh(1), 'Parent'); set(axesh, 'Xlim', [min(D(prev_pos,n-1))-2 max(D(prev_pos,n-1))+2]);
%                 subplot(1,4,2), barh([(D_overlap1(n) + D_overlap2(n) * g(tmpindex)) * 40, beta_SEP(n) * ( tmpindex - optimal_SEP(n) ).^2 / num_lanes,zeros(length(tmpindex),5)],'stacked'); set(gca,'YDir','reverse','XTickLabel',[],'YTickLabel',[]);
%                 subplot(1,4,3), barh(repmat([0,0,2 / num_lanes * f(rawindex,:) * alpha_ideal(n,:)', -PEAK_WEIGHT * peak_bonus(rawindex), -num_lanes * peak_bonus_matrix(rawindex,:) * alpha_ideal(n,:)',0,0] ...
%                     ,length(tmpindex),1),'stacked'); set(gca,'YDir','reverse','XTickLabel',[],'YTickLabel',[]);
%                 subplot(1,4,4), barsh = barh([zeros(length(tmpindex),6), D_test(prev_pos)'],'stacked'); set(gca,'YDir','reverse','XTickLabel',[],'YTickLabel',ttttt);
%                 axesh = get(barsh(1), 'Parent'); set(axesh, 'Xlim', [min(D_test(prev_pos))-2 max(D_test(prev_pos))+2]);
%             end
%         end
    
        [ D(rawindex,n), prev_pos_best(rawindex,n) ] = min( D_test );
        
        D_test(prev_pos) = 10 * BIG_NUMBER;
        
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Backtrack
xsel_fit = [];
[~, xsel_fit(N) ] = min( D(:,N) );
if ( END_POS > 0 ), xsel_fit(N) = END_POS; end;

for n = (N-1) : -1 : 1
    xsel_fit(n) = prev_pos_best( xsel_fit(n+1), n+1 );
end

%figure; image( prev_pos_best )
%figure; image( D/1000 )











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


