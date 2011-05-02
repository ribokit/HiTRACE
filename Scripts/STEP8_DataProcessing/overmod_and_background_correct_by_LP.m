function [ area_correct_bsub, darea_correct_bsub ] = overmod_and_background_correct_by_LP( area_peak, backgd_cols, normbins, area_pred, darea_peak );
% OVERMOD_AND_BACKGROUND_CORRECT_BY_LP: correct for attenuation due to stops in reverse transcription and background correction
% U ses linear programming for effective optimization and (any) guess for what the band intensities 'should' look like (uniform is OK).
%
%  [ area_correct_bsub, darea_correct_bsub ] = overmod_and_background_correct_by_LP( area_peak, backgd_cols, normbins, area_pred, darea_peak );
%
% Inputs:
%  area_peak   = quantitated band intensities for one or more traces. Must include at least one 'background' 
%                     (no-modification) control.
%  backgd_cols = indices corresponding to control(s). If more than one, average is used as background estimate.
%  normbins    = which bands can be used to normalize data. Typically something like [20: size(area_peak,1)-20] 
%  area_pred   = predicted band areas if there was no attenuation or background; usually 0's and 1's. Put all 1's if 
%                  no information -- assuming 'uniform' bands.
%  darea_peak  = error on band intensities [optional]
%
% Output:
%
%  area_correct_bsub  = normalized, background subtracted areas
%  darea_correct_bsub = absolute value of amount background subtracted; OR, if darea_peak defined, appropriately scaled 
%                          errors added in quadrature of signal and background.
%
% (C) R. Das 2010-2011.

goodbins = normbins;
area_peak_in = area_peak;
area_peak = quick_norm( area_peak, normbins );

DAREA_PEAK_DEFINED = 1;
if ~exist( 'darea_peak' ); 
  darea_peak = area_peak * 0; 
  DAREA_PEAK_DEFINED = 0;
end;

for k = 1:size( area_peak, 2 ); darea_peak(:,k) = darea_peak(:,k) * sum( area_peak(:,k) )/sum( area_peak_in(:,k) ); end;

backgd_estimate = mean( area_peak(:, backgd_cols), 2 );
dbackgd_estimate = sum( darea_peak( :, backgd_cols ), 2 );

b = backgd_estimate( goodbins );

% grid search over potential over modification corrections -- probably could be sped up through fminsearch or similar.
overmod_correct = [0.0:0.05:2.0];


for k = 1: size( area_peak,2 )
  
  sum_abs_deviation = [];
  params = [];

  % this inner loop can be totally parallelized!
  if exist( 'matlabpool' )
    if matlabpool( 'size' ) == 0 ;   res = findResource; matlabpool( res.ClusterSize ); end    
    parfor m = 1:length(overmod_correct)
       [area_correct(:,m), params(:,m), sum_abs_deviation(m) ] = overmod_correct_inner_loop( k, m, overmod_correct, goodbins, area_peak, area_pred, backgd_estimate);
    end
  else   
    for m = 1:length(overmod_correct)
      [area_correct(:,m), params(:,m), sum_abs_deviation(m) ] = overmod_correct_inner_loop( k, m, overmod_correct, goodbins, area_peak, area_pred, backgd_estimate);
    end
  end


  [dummy, best_index] = min( sum_abs_deviation );

  area_correct_bsub(:,k) = area_correct(:,best_index) - params(1,best_index)*backgd_estimate;


  if DAREA_PEAK_DEFINED
    darea_correct = ( area_correct(:,best_index) ./ area_peak(:,k) ) .* darea_peak(:,k);
    darea_correct_bsub(:,k) = sqrt( darea_correct.^2 + (params(1,best_index)*dbackgd_estimate).^2 );
  else
    darea_correct_bsub(:,k) = params(1,best_index)* backgd_estimate; % assume error is comparable to amount subtracted here.
  end
  
  %if( params(2,best_index) > 0.0 )
   % area_correct_bsub(:,k) = area_correct_bsub(:,k)/params(2, best_index );
  %elseif( params(1,best_index) > 0.0 )
  %area_correct_bsub(:,k) = area_correct_bsub(:,k)/params(1, best_index );
  %end
  
  fprintf( '%3d. Rate of modification: %8.3f. Background norm:  %8.3f.  Signal strength: %8.3f\n', k, ...
	   overmod_correct(best_index), params(1,best_index), params(2,best_index) );

  plot( overmod_correct, sum_abs_deviation,'k' );
  %pause;
  
end

clf;
subplot(1,3,1)
image( 20*area_peak );

subplot(1,3,2)
image( 20*area_correct_bsub );

subplot(1,3,3)
image( 40*area_pred );

colormap( 1- gray(100));
fprintf('\n');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function    [area_correct, params, sum_abs_deviation ] = overmod_correct_inner_loop( k, m, overmod_correct,goodbins,area_peak,area_pred,backgd_estimate);
area_correct = apply_overmod_correction( area_peak(:,k), overmod_correct(m) );
s = area_correct(goodbins);
y = area_pred(goodbins,k);
b = backgd_estimate( goodbins );
    
[ params, sum_abs_deviation ] = background_fit_LP( s, y, b );
