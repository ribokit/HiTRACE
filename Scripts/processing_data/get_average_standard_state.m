function [d_mean_final, d_err_final] = get_average_standard_state( d_mean, d_err );
% GET_AVERAGE_STANDARD_STATE
%
%  [d_mean_final, d_err_final] = get_average_standard_state( d_mean, d_err );
%
%  combines data from replicate experiments
%  scales up errors if scatter between data sets is higher than estimated error.
%
% INPUTS
%  d_mean = cell of input reactivities
%  d_err  = cell of input reactivity errors.
%
% OUTPUTS
%  d_mean_final = final averaged value 
%  d_err_final  = final comined error
%
if nargin == 0;  help( mfilename ); return; end;

d_mean_final = 0 * d_mean{1};
d_err_final  = 0 * d_mean{1};


% why iterative? 
% need to scale up errors if there is evidence for higher scatter between data sets...
% these are the alpha factors
%alpha = ones( 1, length( d_mean ) );
%NITER = 3;
%for n  = 1 : NITER    
%  for i = 1:length( d_mean )
%    d_mean_final = d_mean_final + d_mean{i}./( alpha(i) * d_err{i}).^2;
%    d_err_final  = d_err_final  + 1./( alpha(i) * d_err{i}).^2;
%  end  
%  d_mean_final = d_mean_final ./ d_err_final;
%  d_err_final = sqrt( 1 ./ d_err_final );
%  for i = 1:length( d_mean )
%    alpha(i) = alpha(i) * mean( ( d_mean_final - d_mean{i} ).^2 ./ d_err{i}.^2 );
%  end
%  alpha = max( alpha, 1.0 );
%end

for i = 1:length( d_mean )
  d_mean_final = d_mean_final + d_mean{i}./d_err{i}.^2;
  d_err_final  = d_err_final  + 1./d_err{i}.^2;
end  
d_mean_final = d_mean_final ./ d_err_final;
d_err_final = sqrt( 1 ./ d_err_final );

% use a sliding window...
window_size = 4;
n = length( d_mean_final );
alpha = 1 + 0*d_mean_final;
for i = 1:n
  window_start = min( max( i-(window_size/2), 1 ),  n - window_size );
  window = window_start + [1:window_size] - 1;
  
  alpha_window = 0;
  for j = 1:length( d_mean )
    alpha_window = alpha_window + mean( (d_mean_final( window ) - d_mean{j}(window) ).^2 ./ d_err{j}(window).^2 );
  end
  alpha_window = sqrt( alpha_window / length( d_mean ) );
  alpha(i) = alpha_window;

end
%plot( alpha );  ylim([0 5]); pause;
alpha = max( alpha, 1.0 );
d_err_final = d_err_final .* alpha;

