function dm = get_weighted_correlation_coefficient( reactivity, reactivity_error );

% this can be optimized later...
N = size( reactivity, 2);

weights = 1./reactivity_error.^2;
reactivity_mean = sum( reactivity.*weights )./sum( weights );
reactivity_shift = reactivity - repmat( reactivity_mean, size(reactivity,1), 1 );

% how about no shifting?
%reactivity_shift = reactivity;% - repmat( reactivity_mean, size(reactivity,1), 1 );

for i = 1:N
  for j = 1:N   
    weights = 1./( reactivity_error(:,i).^2 + reactivity_error(:,j).^2 );
    dm(i,j) = sum( weights .* reactivity_shift(:,i) .* reactivity_shift(:,j) );
    %dm(i,j) = sum( weights .* reactivity_shift(:,i) .* reactivity(:,j) );
    dm(i,j) = dm(i,j) / sqrt( sum( weights .* reactivity_shift(:,i).^2) );
    dm(i,j) = dm(i,j) / sqrt( sum( weights .* reactivity_shift(:,j).^2) );    
  end
end

dm = 1 - dm;
for i = 1:N;
  dm(i,i) = 0;
end
dm = squareform( dm );  % in the right format for MATLAB hierarchical clustering