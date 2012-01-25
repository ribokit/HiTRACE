function [ area_correct_bsub, darea_correct_bsub ] = overmod_and_background_correct_logL( area_peak, backgd_cols, normbins,  area_pred, area_peak_error );
% OVERMOD_AND_BACKGROUND_CORRECT_LOGL
%
%  [ area_correct_bsub, darea_correct_bsub ] = overmod_and_background_correct_by_LP( area_peak, backgd_cols );
%
% Inputs:
%  area_peak   = quantitated band intensities for one or more traces. Must include at least one 'background' 
%                     (no-modification) control.
%  backgd_cols = indices corresponding to control(s). If more than one, average is used as background estimate.
%  normbins    = [ OPTIONAL ] which bands can be used to normalize data. Typically something like [20: size(area_peak,1)-20]   
%  area_pred   = [ OPTIONAL ] predicted band areas if there was no attenuation or background; usually 0's and 1's. Put all 1's if 
%                  no information -- assuming 'uniform' bands.
%
% Output:
%
%  area_correct_bsub  = normalized, background subtracted areas
%  darea_correct_bsub = absolute value of amount background subtracted; OR, if darea_peak defined, appropriately scaled 
%                          errors added in quadrature of signal and background.
%
% (C) R. Das 2011.


if ~exist( 'normbins' )
  normbins = [1:size( area_peak,1) ];
end
if ~exist( 'area_pred' )
  area_pred = [];
end
if ~exist( 'area_peak_error' )
  area_peak_error = abs(area_peak) * 0.2;
end

b = mean( area_peak( :, backgd_cols ), 2 );
db = sqrt( sum( area_peak_error(:,backgd_cols ).^2, 2 )/length( backgd_cols) );

for k = 1: size( area_peak,2 )

  %if sum( area_pred(:,k) ) == 0  % background column
  %  if ( length(backgd_cols) == 1 & k == backgd_cols )
  %    area_correct_bsub(:,k) = 0 * b;   alpha = 0.0; beta = 1.0; % solution is degenerate
  %  else
  %    [ area_correct_bsub(:,k), alpha, beta ] = backsub_and_norm_logL( area_peak(:,k), b, normbins, 1/0.05, 1/0.05);
  %  end
  %else
  %  [ area_correct_bsub(:,k), alpha, beta ] = overmod_wrapper_logL( area_peak(:,k), b, normbins );
  %end

  area_pred_to_use = [];
  if ~isempty( area_pred ); area_pred_to_use = area_pred(:,k); end;
  if ~isempty( find( k == backgd_cols ) );  area_pred_to_use = 0 * area_peak(:,k); end; 
  [ area_correct_bsub(:,k), alpha, beta ] = overmod_wrapper_logL( area_peak(:,k), b, normbins, area_pred_to_use );
  darea_correct_bsub(:,k) = alpha * sqrt( area_peak_error(:,k).^2 + ( beta * db).^2);

  %plot( overmod_correct, sum_abs_deviation,'k' );
  %pause;
  
end