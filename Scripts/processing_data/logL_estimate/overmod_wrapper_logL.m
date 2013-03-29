function [ s_bsub,  alpha, beta, L, s_bsub_all, overmod ]  = overmod_wrapper_logL(  s, b, normbins, area_pred, overmod_specified, PLOT_STUFF );
% [ s_bsub,  alpha, beta, L, s_bsub_all, overmod ]  = overmod_wrapper_logL(  s, b, normbins, area_pred, overmod_specified, PLOT_STUFF );
%

if nargin == 0;  help( mfilename ); return; end;

if ~exist( 'normbins' ) | isempty( normbins); normbins = [1:length(b)]; end;
if ~exist( 'area_pred' ); area_pred = 0*b - 1; end;
if ~exist( 'overmod_specified' ); overmod_specified = -1; end;
if ~exist( 'PLOT_STUFF' ); PLOT_STUFF = 0; end;

overmod = [0.0:0.1:4.0];
%overmod = [0.0:0.05:2.0];
if ( overmod_specified >= 0 )
  overmod = [overmod_specified];
end

n = length( overmod );

if abs( sum( s - b ) ) == 0.0  % signal = background?
  s_bsub = s - b;
  alpha = 1.0;
  beta = 1.0;
  L = 0.0;
  s_bsub_all = [];
  fprintf( ' Rate of modification: %8.3f. Background norm:  %8.3f.  Signal strength: %8.3f\n', ...
	   0.0, beta, alpha )
  return;
end


s_save = s;
b_save = b;
s = s( normbins );
b = b( normbins );
if ~isempty( area_pred ); area_pred = area_pred( normbins ); end;

s_bsub_all = zeros( length( s ), n );
alpha = zeros(1,n);
beta = zeros(1,n);
L = zeros(1,n); % this is -log Likelihood

% this inner loop can be totally parallelized!
if parallelization_exists()
  if matlabpool( 'size' ) == 0 ;   res = findResource; matlabpool( res.ClusterSize ); end    
  parfor i = 1:n
    [s_correct, correction ] = apply_overmod_correction_EXACT( s, overmod(i) ) ;
    %[s_bsub_all(:,i), alpha(i), beta(i), L(i) ] =  backsub_and_norm_logL_with_cutoffs( s_correct, b );
    [s_bsub_all(:,i), alpha(i), beta(i), L(i) ] =  backsub_and_norm_logL( s_correct, b, [], area_pred );
    L(i) = L(i) - sum( log( correction ) );
  end
else
  for i = 1:n
    [ s_correct, correction ] = apply_overmod_correction_EXACT( s, overmod(i) ) ;
    %[s_bsub_all(:,i), alpha(i), beta(i), L(i) ] =  backsub_and_norm_logL_with_cutoffs( s, b );
    [s_bsub_all(:,i), alpha(i), beta(i), L(i) ] =  backsub_and_norm_logL( s_correct, b, [], area_pred );
    L(i) = L(i) - sum( log( correction ) );
  end
end

[ bestL, best_idx ] = min( L );

alpha = alpha( best_idx );
beta  = beta( best_idx );

s_correct = apply_overmod_correction_EXACT( s_save, overmod(best_idx) );
s_bsub = alpha * ( s_correct - beta * b_save );

%[overmod( best_idx ) alpha( best_idx)  beta( best_idx) ]

fprintf( ' Rate of modification: %8.3f. Background norm:  %8.3f.  Signal strength: %8.3f\n', ...
	 overmod(best_idx), beta, alpha )

if PLOT_STUFF
  subplot(2,1,1)
  plot( overmod, L );
  
  subplot(2,1,2)
  plot( [alpha*s_correct, alpha*beta*b_save, s_bsub] )
end
