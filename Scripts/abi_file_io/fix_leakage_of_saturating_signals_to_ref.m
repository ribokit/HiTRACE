function da_fix = fix_leakage_of_saturating_signals_to_ref(d, da)
%
%
%

da_fix = da;

for i = 1:size(d,2);
  max_d = max( d(:,i) );
  % this is kind of extreme...
  cutoff =  0.90 * max_d; 
  
  badpoints = find( d(:,i) > cutoff );
  da_fix( badpoints, i ) = 0.0;
end
