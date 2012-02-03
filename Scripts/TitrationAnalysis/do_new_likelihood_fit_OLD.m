function  [ logL, pred_fit, lane_normalization, sigma_at_each_residue, C_state ] = ...
    do_new_likelihood_fit( data, conc, K1, n, fit_type, lane_normalization, C_state_in, beta_C  );

if ~exist( 'beta_C' ); beta_C = 0.1; end;
if ~exist( 'fit_type' ); fit_type = 'hill'; end;

%Expect a minimum of ~10% error due to systematics.
SIGMIN_FRAC = 0.1; 

%pred = get_fraction_extended_Hill( conc, K1, n );
%f = [ 1-pred; pred ]; %fraction in state.

f = feval( fit_type, conc, K1, n );

numres = size( data, 1 );
numconc = size( data, 2);

numiter = 4;
sigma_normalization = 1.0;

if ~exist( 'lane_normalization' ) | isempty( lane_normalization)
  lane_normalization = ones(1, numconc); 
  numiter = 1; % iteration is not working??? 
end;
if ~exist( 'C_state_in' )
  C_state_in = [];
end
sigma_at_each_residue = ones( 1, numres );

for n = 1:numiter

  [pred_fit, pred_fit_before_scaling, C_state ] = do_linear_fit_vs_conc( data, f, lane_normalization, C_state_in, sigma_at_each_residue, beta_C );

  sigma_at_each_residue = get_sigma_at_each_residue( data, pred_fit, SIGMIN_FRAC );
  if ( n < numiter )
    lane_normalization = normalize_lanes( data, pred_fit, pred_fit_before_scaling, sigma_normalization, sigma_at_each_residue );
  end

  sigma_normalization = sqrt( mean(( lane_normalization - 1 ).^2) );

  %sigma_normalization =std( lane_normalization ); 
  
  %logL = - numconc * sum( log( sigma_at_each_residue ) ) ...
  %- numconc * log( sigma_normalization );

  %plot( lane_normalization ); hold on
  %pause;
  
end

logL = - numconc * sum( log( sigma_at_each_residue ) );

if ( sigma_normalization > 0 )  logL = logL - numconc * numres * log( sigma_normalization ); end;

if exist( 'C_in' ) logL = logL - beta_C * sum( sum( (( C_in - C_state ) * diag( 1./sigma_normalization)).^2 ) ); end;

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pred_fit, pred_fit_before_scaling, C ] = do_linear_fit_vs_conc( data, f, lane_normalization, C_in, sigma_at_each_residue, beta );

numres = size( data, 1 );
numconc = size( data, 2);
num_states = size( f, 1 );
[lane_norm_grid, dummy] = meshgrid( lane_normalization, 1:numres );

JUST_SCALE_CIN = 0;
if exist( 'C_in' ) & ~isempty( C_in ) & JUST_SCALE_CIN

  % this is a special case -- we 'know' what the states look like, and just
   % need to scale them proportionally.
  pred_fit = lane_norm_grid .* ( C_in' * f );
  for a = 1:num_states
    pred_fit_contribution = lane_norm_grid .* (C_in( a, : )' * f(a, :));
    %clf; plot( pred_fit ); hold on; plot( pred_fit_contribution, 'r'); pause;
    kappa(a) = sum(sum(data .* pred_fit_contribution)) / sum(sum(pred_fit .* pred_fit_contribution));
  end
  C = diag(kappa) * C_in;

else
  A = f * diag( lane_normalization.^2 ) * f';
  % reactivity of each state.
  B =  data*diag(lane_normalization)*f';

  if exist( 'C_in' ) & ~isempty( C_in )
    A = A + beta;
    B = B + beta * C_in';
  end
  
  C = A\B';
  
  % to enforce that C is positive...
  %for m = 1:size( B, 1 );
  %  C(:,m) = lsqnonneg( A, B(m,:)' );
  %end
end

pred_fit_before_scaling = C' * f;
pred_fit = lane_norm_grid .* pred_fit_before_scaling; 

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sigma_at_each_residue = get_sigma_at_each_residue( data, pred_fit, SIGMIN_FRAC );

deviation_fit = (pred_fit - data);
sigma_at_each_residue2 = mean( deviation_fit.^2, 2 ) + ...
    (SIGMIN_FRAC *mean( data,2)).^2   ;
sigma_at_each_residue = sqrt( sigma_at_each_residue2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ lane_normalization, sigma_at_each_residue] = ...
    normalize_lanes( data, pred_fit, pred_fit_before_scaling, ...
		     sigma_normalization, sigma_at_each_residue );

numres = size( data, 1);
numconc = size( data, 2);

deviation_fit = (pred_fit - data);
covariance_inverse = diag( 1./sigma_at_each_residue.^2);

USE_SOPHISTICATED_COVARIANCE = 0;
if (USE_SOPHISTICATED_COVARIANCE)
  covariance_matrix = (deviation_fit * deviation_fit')/ ( numconc );
  covariance_matrix_apodized = apodize(covariance_matrix);
  %covariance_matrix_apodized = covariance_matrix + eye( numres ) * SIGMIN^2;
  covariance_inverse = inv( covariance_matrix_apodized );
end

%sigma_alpha = 0.1;
%reg_factor = numres * (1/(sigma_alpha^2));

%reg_factor = 1/(sigma_normalization^2); % this is wrong
reg_factor = numres/(sigma_normalization^2); % this is right.

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
