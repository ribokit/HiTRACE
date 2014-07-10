function [final_reactivity, final_error, flags ] = average_data_filter_outliers( reactivities, guessed_errors, seqpos, disallow_outlier_res ); 
% AVERAGE_DATA_FILTER_OUTLIERS
% [final_reactivity, final_error, flags ] = average_data_filter_outliers( reactivities, guessed_errors, seqpos, disallow_outlier_res ); 
%
% Takes traces and initial error estimates; figures out outlier points and even outlier traces;
%  and then returns reasonable final values and error estimates.
%
% Inputs:
%  reactivities   = Matrix of input traces
%  guessed_errors = Matrix of corresponding error estimates. These are used to weight the average and
%                     filter outliers, but final_error depends on the actual scatter seen in (non-outlier) points.
%  seqpos         = [Optional] conventional sequence numbering -- for warnings.
%  disallow_outlier_res = [Optional] define residues to be excluded from outlier filtering 
%
% Outputs:
%  final_reactivity = final averaged reactivity, weighted by 1/error^2.
%  final_error      = error on final reactivity, from standard deviation across traces / sqrt(N).
%  flags            = matrix of 0 and 1's, with 1 flagging problem points.
%
% (C) R. Das, T. Mann, Stanford University, 2013.

POINT_ERROR_RATIO_CUTOFF = 5.0;
PROFILE_ERROR_RATIO_CUTOFF = 2.5;

if ( nargin < 2 ) help( mfilename ); return; end;

N = size( reactivities, 2 );
L = size( reactivities, 1 );
if ~exist( 'seqpos','var' ); seqpos = [1:L]; end;
if ~exist( 'disallow_outlier_res','var' ); disallow_outlier_res = []; end;
disallow_outlier_pos = disallow_outlier_res - min(seqpos) + 1;

flags = ones(L,N);
final_err = zeros(N,1);
final_reactivity = zeros(N,1);

NITER = 5;
for i = 1:L
  gp = 1:N;

  weights = max( 1./guessed_errors(i,:), 0 ) ;

  for n = 1:NITER
    m = sum( reactivities( i, gp ) .*weights(gp) ) / sum( weights(gp) ) ;   % calculate weighted sum
    dev = abs(reactivities(i,:) - m );  % 
    if any( i == disallow_outlier_pos );
        fprintf( 'Not excluding residue %d\n', seqpos(i) );
    else
        gp = find( dev < POINT_ERROR_RATIO_CUTOFF * guessed_errors(i,:)  );
    end
  end  
  flags(i,gp) = 0;
  
  if length( gp ) == 0; % whoa that's a big problem.
    fprintf( 'Residue %d has way more scatter than input guessed_errors!\n', seqpos(i) );
    gp = [1:N];
  end
  
  m = sum( reactivities( i, gp ) .* weights(gp) ) / sum( weights(gp) );
  final_reactivity(i) = m; 
  
  if length( gp ) < 2; gp = [1:N]; end;
  s = std( reactivities(i, gp) );
  final_error(i) = max( s/sqrt( length(gp)), 0 );
end
final_error = final_error';

% now look over all input reactivity profiles -- are there some that are just crazy?
bad_traces = [];
for j = 1:N
  mean_rel_error(j) = mean( max(abs(reactivities(:,j) - final_reactivity) ./ final_error,0) );
end
for j = 1:N  
  %[mean_rel_error(j) median( mean_rel_error )]
  if ( mean_rel_error(j) > (PROFILE_ERROR_RATIO_CUTOFF * median( mean_rel_error )) )
    bad_traces = [bad_traces, j];
  end
end
good_traces = setdiff( [1:N], bad_traces );

if length( good_traces ) > 3 & length( bad_traces ) > 0
  fprintf( ['\nFollowing traces look way off -- will not use them for averaging: ',num2str(bad_traces),'\n\n'] );
  flags(:,bad_traces ) = 1.0;

  good_traces = setdiff( [1:N], bad_traces );
  [final_reactivity, final_error, flags(:,good_traces)] = average_data_filter_outliers( reactivities(:,good_traces), guessed_errors(:,good_traces), seqpos, disallow_outlier_res );
end