function dm = get_chi_squared_dm( reactivity, reactivity_error );

% this can be optimized later...
N = size( reactivity, 2);
for i = 1:N
  for j = 1:N
    chi2(i,j) = sum( (reactivity(:,i) - reactivity(:,j)).^2 ./ ( reactivity_error(:,i).^2 + reactivity_error(:,j).^2) );
  end
end

dm = sqrt( chi2 ) / N;
dm = squareform( dm );  % in the right format for MATLAB hierarchical clustering