function [d_collate, bounds] = collate_data( data, CONCAT_TRACES );

if ~exist( 'CONCAT_TRACES', 'var' ) CONCAT_TRACES = 1; end;
bounds= [];
d_collate = [];
for i = 1:length( data )
  if CONCAT_TRACES
    d_collate = [ d_collate; data{i} ];
    bounds(i) = size( d_collate, 1 );
  else
    d_collate = [ d_collate, data{i} ];
    bounds(i) = size( d_collate, 2 );
  end
end