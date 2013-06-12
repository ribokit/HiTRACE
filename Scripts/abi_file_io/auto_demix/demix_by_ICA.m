function data_out = demix_by_ICA( data, num_components );
% DEMIX_BY_ICA: Totally experimental scripts for automatically 
%     carrying out  leakage correction through Independent Component Analysis.
%
%  data_out = demix_by_ICA( data, num_components );
%
% data = cell of matrices with traces, as would be output by plot_ABI_runs. 
%         Each matrix has, e.g., four 
%         columns for the four color channels.
%
% num_components = number of actual signal components (default 2).
%
% This doesn't really work unless input data has been subjected to
% baseline_subtract_smooth. Even then, some weirdness happens due to
% saturated peaks -- need to work on that.
% 
% (C) R. Das, 2013
%

if ~exist( 'num_components' ) num_components = 2; end;
[d_orig, bounds_side_by_side_orig] = collate_data( data, 0 );
[d, bounds] = collate_data( data, 1 );

dx = fastica( d', 'approach','symm', 'lastEig', num_components ); plot( dx' )
dx = dx';

% flip signs!
dx = flip_signs( dx );

% figure out reordering that best matches original data.
for i = 1:size( d, 2 );
  for j = 1:size( dx, 2 );
    c = corrcoef( d(:,i), dx(:,j) );
    m(i,j) = c(1,2);
  end
end

correspondence = [];
for j = 1:size( dx, 2 ) % reseparate.
  [~,sortidx] = sort( m(:,j) );
  i = sortidx(end);
  while ~all( find( i ~= correspondence ) )
    sortidx = sortidx( 1: end-1 ); 
    i = sortidx( end );
  end
  correspondence(j) = i;
end

for j = 1:size( dx, 2 ) % reseparate.
  dx(:,j) = dx(:,j) * mean( d(:, correspondence(j) ) )/ mean( dx(:,j) );
end


[~,reorder] = sort( correspondence );
dx = dx(:,reorder);



data_out = separate_data( dx, bounds );

[d_out, bounds_side_by_side_out] = collate_data( data_out, 0 );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clf
make_colormap;
scalefactor = 20/mean(mean( d_orig))
subplot(1,2,1);
image( d_orig*scalefactor )
make_lines( bounds_side_by_side_orig,'b',1 );
set(gca,'tickdir','out');
title( 'input' );

subplot(1,2,2);
image( d_out*scalefactor )
make_lines( bounds_side_by_side_out,'b',1 );
set(gca,'tickdir','out');
title( 'after ICA' );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dx = flip_signs( dx )
for i = 1:size( dx, 2 );
  dx(:,i) = dx(:,i) * sign( mean(dx(:,i)) );
end
