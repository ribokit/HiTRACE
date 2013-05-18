function I_norm = boxcart_normalize( I_data, window_size );

if ~exist( 'window_size' ), window_size = 50; end;
I_norm = 0 * I_data;
n = size( I_data, 1 );

min_norm_val = mean( std( I_data ) )/2;

for i = 1:n
  idx_min = max(1,   round(i - window_size/2)  );
  idx_max = min(n,   round(i + window_size/2)  );
  mean_val = mean( mean( I_data( idx_min:idx_max, : ) ) );
  mean_val = max( min_norm_val, mean_val );
  I_norm(i,:) = I_data(i,:)/mean_val;
end
