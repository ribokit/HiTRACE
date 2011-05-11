function  [ logL, pred_fit, lane_normalization, sigma_at_each_residue,...
	    data_rescale, pred ] = ...
    do_new_likelihood_fit( data, conc, resnum, K1, n );

%Expect a minimum of ~10% error due to systematics.
SIGMIN_FRAC = 0.1;

pred = get_fraction_extended_Hill( conc, K1, n );

numres = size( data, 1 );
numconc = size( data, 2);

numiter = 2;
sigma_normalization = 1.0;
lane_normalization = ones(1, numconc);


for n = 1:numiter
  
  [pred_fit, pred_fit_before_scaling, data_rescale ] = do_linear_fit_vs_conc( data, pred, lane_normalization );

  [lane_normalization, sigma_at_each_residue] = ...
      normalize_lanes( data, pred_fit, pred_fit_before_scaling, ...
		       sigma_normalization, SIGMIN_FRAC );

  sigma_normalization = sqrt( mean(( lane_normalization - 1 ).^2) ); 
  
%More correct:
  logL = - numconc * sum( log( sigma_at_each_residue ) ) ...
	 - numconc * log( sigma_normalization );

end

%Hmmm
%logL = get_likelihood( data, pred_fit, SIGMIN );
%logL = logL - numconc * log( sigma_normalization );

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pred = get_fraction_extended_Hill( conc, K1, n );

pred = (conc/K1).^n ./ (1 + (conc/K1).^n );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pred_fit, pred_fit_before_scaling, data_rescale] = do_linear_fit_vs_conc( data, pred, lane_normalization );

numres = size( data, 1 );
numconc = size( data, 2);


if ~exist( 'lane_normalization' )
  lane_normalization = ones( 1, numconc );
end

%This is in crazy matlab-ese to make it fast...
[lane_norm_grid, dummy] = meshgrid( lane_normalization, 1:numres );
data_unscale = data./lane_norm_grid;

scaled_lane_normalization2 = lane_normalization.^2 / mean( lane_normalization.^2 );						 
[dummy, data_unscale_avg_grid] = meshgrid( 1:numconc, ...
					  (data_unscale*scaled_lane_normalization2')/numconc );

data_unscale_center = data_unscale - data_unscale_avg_grid;


pred_center = pred - mean( pred.*scaled_lane_normalization2 );
[pred_center_grid, dummy] = meshgrid( pred_center, 1:numres );
[scaled_lane_norm2_grid, dummy] = meshgrid( scaled_lane_normalization2, 1:numres );

rescale_factor = sum(data_unscale_center.*pred_center_grid.*scaled_lane_norm2_grid, 2)./...
    sum(pred_center_grid.*pred_center_grid.*scaled_lane_norm2_grid, 2);

[dummy, rescale_factor_grid] = meshgrid( 1:numconc, rescale_factor );

pred_fit_before_scaling = ( pred_center_grid .* rescale_factor_grid ...
			    + data_unscale_avg_grid ) ;

pred_fit = lane_norm_grid .* pred_fit_before_scaling; 

%useful for plotting:
data_rescale = ( data_unscale_center ./ rescale_factor_grid ) + mean( pred.*scaled_lane_normalization2 );

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ lane_normalization, sigma_at_each_residue] = ...
    normalize_lanes( data, pred_fit, pred_fit_before_scaling, ...
		     sigma_normalization, SIGMIN_FRAC );

numres = size( data, 1);
numconc = size( data, 2);

deviation_fit = (pred_fit - data);
sigma_at_each_residue2 = mean( deviation_fit.^2, 2 ) + ...
    (SIGMIN_FRAC *mean( data,2)).^2   ;
covariance_inverse = diag( 1./sigma_at_each_residue2);

USE_SOPHISTICATED_COVARIANCE = 0;
if (USE_SOPHISTICATED_COVARIANCE)
  covariance_matrix = (deviation_fit * deviation_fit')/ ( numconc );
  covariance_matrix_apodized = apodize(covariance_matrix);
  %covariance_matrix_apodized = covariance_matrix + eye( numres ) * SIGMIN^2;
  covariance_inverse = inv( covariance_matrix_apodized );
end

reg_factor = 1/(sigma_normalization^2);

for j = 1:numconc
  
%  lane_normalization( j ) = ...
%      ( sum(data(:,j).*pred_fit(:,j)./ ...
%	  (sigma_at_each_residue.^2) )  + reg_factor )/ ...
%      ( sum(pred_fit(:,j).*pred_fit(:,j)./ ...
%	  (sigma_at_each_residue.^2) ) + reg_factor);
    lane_normalization( j ) = ...
      ( data(:,j)'*covariance_inverse*pred_fit_before_scaling(:,j) + reg_factor )/ ...
      ( pred_fit_before_scaling(:,j)'*covariance_inverse*pred_fit_before_scaling(:,j) + reg_factor);

end

sigma_at_each_residue = sqrt( sigma_at_each_residue2);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function logL = get_likelihood( data, pred_fit, SIGMIN );

numres = size( data, 1 );
numconc = size( data, 2);

logL = 0.0;
if find( isnan( pred_fit  ) ) 
  return
end

deviation_fit = (pred_fit - data);

covariance_matrix = (deviation_fit * deviation_fit') / size( data, 1);

covariance_matrix_apodized = apodize( covariance_matrix );
%covariance_matrix_apodized = covariance_matrix + SIGMIN^2 * eye(numres);

[u,s,v] = svd( covariance_matrix_apodized, 0 );
logL = -0.5 * numconc * sum( log( diag(s) + SIGMIN^2 ) );
%logL = -log(det( covariance_matrix ));

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function covariance_matrix_apodized = apodize( covariance_matrix );
numres = size( covariance_matrix, 1);
[i_grid,j_grid] = meshgrid(1:numres, 1:numres);

MAX_CORRELATION = 1000.0;
covariance_matrix_apodized = covariance_matrix .* ...
    exp( - abs(i_grid-j_grid) / MAX_CORRELATION );
