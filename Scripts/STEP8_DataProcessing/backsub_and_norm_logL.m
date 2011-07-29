function  [x, alpha, beta, L ] = backsub_and_norm_logL( s, b, normbins, area_pred );
% [x, alpha, beta, L ] = backsub_and_norm_logL( s, b, normbins, area_pred );
%
%  Likelihood-based estimation of background/normalization factor.
%  Corrected data is:
%
%     x = alpha * ( s - beta * b )
%
%  Likelihood function takes the form of  
%
%    alpha * exp( -F_plus  * x ) (for x > x0)
%    alpha * exp( -F_minus * x ) (for x < x0)
%
% where F_plus = 5 [a rough exponential distribution], F_minus = 25 [strong penalties for negative values]
%  and x0 = 0.0. 
%
%
%  INPUTS
%   s        = measured signal
%   b        = estimated background
%   normbins = [optional] positions over which normalization & background scalings are optimized
%   area_pred= [optional] guesses of which positions are reactive or protected (given as 1 and 0, respectively)
%
%  OUTPUT
%   x        = normalized, background-subtracted data.  x = alpha * ( s - beta * b )
%   alpha    = maximum likelihood estimate of alpha (overall normalization)
%   beta     = maximum likelihood estimate of beta  (background scaling)
%   L        = minus log-likelihood value.
%

PLOT_STUFF = 0;

N = length( s );
alpha = 0.25/mean(s);
beta = 0.3;%(1/alpha) * 0.125/mean(b);

if ~exist( 'normbins' ) | isempty(normbins)
  normbins = [1:N];
end

x0 = zeros( N, 1);  % 'most likely value'
F_plus  = zeros( N, 1); % attenuation factor for signal over the most likely value
F_minus = zeros( N, 1); % attenuation factor for signal under the most likely value

if length( area_pred ) > 0
  for k = 1:N
    if area_pred(k)
      F_plus(k) = 1/0.39;
      F_minus(k) = 1/0.04;
      x0(k) = 0.06;
      %x0(k) = 0.0;
    else
      %F_plus(k) = 1/0.10;
      F_plus(k) = 1/0.04;
      F_minus(k) = 1/0.04;
      x0(k) = 0.0;
    end
  end
else
  F_plus  = 1/0.2  * ones( N, 1);
  F_minus = 1/0.04 * ones( N, 1);
  x0 = 0.06 * ones( N, 1);
  %x0 = 0.0 * ones( N, 1);
end

s_save = s;
b_save = b;
s = s( normbins );
b = b( normbins );
F_plus  = F_plus( normbins );
F_minus = F_minus( normbins );

niter = 10;
for k = 1:niter
  x = alpha* ( s - beta * b );


  xp = find( x > x0 );
  xm = find( x < x0 );


  L = sum( F_plus(xp) .* ( x(xp) - x0(xp) ) )...
      - sum( F_minus(xm) .*  ( x(xm) - x0(xm) ) )...
     - N * log( alpha );    

  if PLOT_STUFF
    [alpha, beta, L]
    plot( [alpha*s, alpha * beta*b, x ] );
    hold on
    plot( xp, x( xp ), '+' );
    plot( xm, x( xm ), 'v' );
    hold off
    pause;      
  end
  
  %alpha = N / ( F_plus  * sum(s( xp )) - F_minus * sum(s( xm )) );
  s_sub = s - beta * b;
  alpha = N / ( sum( F_plus(xp)  .* s_sub( xp )) - sum( F_minus(xm) .* s_sub( xm )) );
    
  % old style.
  vals = (s - x0/alpha)./b;
  [possible_beta, sort_index ] = sort( vals );

  F_balance = [];
  for i = 1:length(vals)
    idx_plus  = sort_index( (i+1): length(vals) );
    idx_minus = sort_index( 1:i );

    b_plus(i)  = sum( F_plus( idx_plus ) .* b( idx_plus)  );
    b_minus(i) = sum( F_minus( idx_minus) .* b( idx_minus) );
    F_balance(i) = b_plus(i) - b_minus(i);
  end
  beta = interp1( F_balance, possible_beta, 0 );

  %beta  = max( beta, 0.0 );
  
  %plot( possible_beta, [ F_balance; F_balance_new] );
  %pause;

  %beta  = N / ( -F_plus  * sum(b( xp )) + F_minus * sum(b( xm )) );
  
  %alpha = max( alpha, 0 );
  beta  = max( beta, 0 );
  
end


x = alpha* ( s_save - beta * b_save );
