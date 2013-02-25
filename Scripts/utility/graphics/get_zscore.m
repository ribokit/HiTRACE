function c = get_zscore( a, ref_cols )

if ~exist( 'ref_cols' ); ref_cols = [1:size(a,2)]; end;

c = [];
a_norm = a;

for j = 1:size(a,1 )
  error = std( a_norm(j, : ) );
  c(j,:) = ( a_norm(j,:) - mean( a_norm(j,ref_cols)) ) /error;
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:size(c,2 )
  c(:,i) = c(:,i) /  std( c(:,i ) );
end

