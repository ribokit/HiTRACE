function I_norm = window_normalize( I_data, window_size );
% WINDOW_NORMALIZE = use moving-window ('window') to normalize data
%
% I_norm = window_normalize( I_data, window_size );
%
% Note that mean is computed ACROSS all traces.
%
% I_data       = input traces
% window_size  = size of window to determine mean [default 50].
%
% (C) R. Das, 2013

if ~exist( 'window_size' ), window_size = 50; end;
I_norm = 0 * I_data;
nres   = size( I_data, 1 );
nlanes = size( I_data, 2 );

% need to remove outliers when computing standard deviation.
NITER = 3;
I_data_flat = reshape( I_data, 1, nres*nlanes);
sigma_cutoff = 5; 
for k = 1:NITER
  s = std( I_data_flat );
  m = mean( I_data_flat );
  gp = find( abs(I_data_flat - m) < sigma_cutoff * s );
  I_data_flat = I_data_flat( gp );
end
min_norm_val = s;

for i = 1:nres
  idx_min = max(1,    round(i - window_size/2)  );
  idx_max = min(nres, round(i + window_size/2)  );
  mean_val = mean( mean( I_data( idx_min:idx_max, : ) ) );
  mean_val = max( min_norm_val, mean_val );
  I_norm(i,:) = I_data(i,:)/mean_val;
end
