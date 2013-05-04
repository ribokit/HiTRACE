function d_offset = remove_offset( d );
% REMOVE_OFFSET
%
% d_offset = remove_offset( d );
%
% Remove constant offset from data. Still under testing.
%

d_offset = d;

for i = 1:size( d, 2 );
  d_offset(:,i) = remove_offset_one_trace( d(:,i) );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d_offset = remove_offset_one_trace( d );

NITER = 10;
d_iter = d;

for n = 1:NITER
  stdev = std( d_iter );
  gp = find( abs( d_iter - mean( d_iter ) ) <  1.5*stdev );
  if length( gp ) > 10
    d_iter = d_iter( gp );
  end
end
d_offset = d - mean( d_iter );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% old - remove this!
function d_offset = remove_offset_one_trace_OLD( d );

stdev = std( d );
fineness = stdev/5;

d_discrete = fineness * round( smooth(d) / fineness );
offset = mode( d_discrete );

d_offset = d - offset;