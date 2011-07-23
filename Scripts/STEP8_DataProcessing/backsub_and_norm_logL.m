function  [g, alpha, beta, L ] = backsub_and_norm_logL( s, b, normbins, area_pred );

PLOT_STUFF = 0;

N = length( s );
alpha = 0.25/mean(s);
beta = 0.3;%(1/alpha) * 0.125/mean(b);

if ~exist( 'normbins' ) | isempty(normbins)
  normbins = [1:N];
end

g0 = zeros( N, 1);  % 'most likely value'
F_plus  = zeros( N, 1); % attenuation factor for signal over the most likely value
F_minus = zeros( N, 1); % attenuation factor for signal under the most likely value

if length( area_pred ) > 0
  for k = 1:N
    if area_pred(k)
      F_plus(k) = 1/0.39;
      F_minus(k) = 1/0.04;
      g0(k) = 0.06;
      %g0(k) = 0.0;
    else
      %F_plus(k) = 1/0.10;
      F_plus(k) = 1/0.04;
      F_minus(k) = 1/0.04;
      g0(k) = 0.0;
    end
  end
else
  F_plus  = 1/0.2  * ones( N, 1);
  F_minus = 1/0.04 * ones( N, 1);
  g0 = 0.06 * ones( N, 1);
  %g0 = 0.0 * ones( N, 1);
end

s_save = s;
b_save = b;
s = s( normbins );
b = b( normbins );
F_plus  = F_plus( normbins );
F_minus = F_minus( normbins );

niter = 10;
for k = 1:niter
  g = alpha* ( s - beta * b );


  gp = find( g > g0 );
  gm = find( g < g0 );


  L = sum( F_plus(gp) .* ( g(gp) - g0(gp) ) )...
      - sum( F_minus(gm) .*  ( g(gm) - g0(gm) ) )...
     - N * log( alpha );    

  if PLOT_STUFF
    [alpha, beta, L]
    plot( [alpha*s, alpha * beta*b, g ] );
    hold on
    plot( gp, g( gp ), '+' );
    plot( gm, g( gm ), 'v' );
    hold off
    pause;      
  end
  
  %alpha = N / ( F_plus  * sum(s( gp )) - F_minus * sum(s( gm )) );
  s_sub = s - beta * b;
  alpha = N / ( sum( F_plus(gp)  .* s_sub( gp )) - sum( F_minus(gm) .* s_sub( gm )) );
    
  % old style.
  vals = (s - g0/alpha)./b;
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

  %beta  = N / ( -F_plus  * sum(b( gp )) + F_minus * sum(b( gm )) );
  
  %alpha = max( alpha, 0 );
  beta  = max( beta, 0 );
  
end


g = alpha* ( s_save - beta * b_save );
