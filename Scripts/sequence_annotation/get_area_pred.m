function area_pred = get_area_pred( sequence, data_types, offset, structure );
% GET_AREA_PRED
%
% [ marks, area_pred, mutpos] = ( sequence, data_types, offset, structure );
%
% 
% Returns information on where bands should show up, based on input 
%  sequence, structure, and what modification/ladder reaction is in each lane.
%
% Used by ANNOTATE_SEQUENCE
%
% sequence   = full sequence
% data_types = cell of tags of modification reactions in each trace, e.g., 
%               {'SHAPE','SHAPE','ddTTP'}
% offset     = value that is added to sequence index to achieve 'historical'/favorite numbering. [default: 0]
% structure  = structure in dot/bracket notation [give as '' if unknown]
%
% OUTPUT:
% area_pred = matrix of zeros and ones corresponding to where bands will 
%              show up.
%
% (C) R. Das, 2013.
%

if nargin == 0;  help( mfilename ); return; end;

seqpos = 1:length(sequence);
area_pred = zeros( length(sequence), length(data_types) );
if ~exist( 'structure' ) | length( structure ) == 0;   for i = 1:length( sequence ); structure = [structure, '.']; end; end
if length( structure ) ~= length( sequence ); fprintf( 'Sequence length must equal structure length!\n'); return; end;


for k = 1:length( data_types )

  switch data_types{k}
   case 'SHAPE'
    for i = 1:length( structure )
      if structure(i) == '.'
	area_pred(i,k) = 1.0;
      end
    end
   case 'DMS'
    for i = 1:length( structure )
      if structure(i) == '.' & ( sequence(i) == 'A'  | sequence(i) == 'C'  ) 
	area_pred(i,k) = 1.0;
      end
    end
   case 'CMCT'
    for i = 1:length( structure )
      if structure(i) == '.' & ( sequence(i) == 'G'  | sequence(i) == 'U'  | sequence(i) == 'T'  ) 
	area_pred(i,k) = 1.0;
      end
    end
   case 'UV'
    for i = 2:length( structure )
      if ( sequence(i) == 'C'  | sequence(i) == 'U'  ) & ...
	    ( sequence(i-1) == 'C'  | sequence(i-1) == 'U'  ) 
	area_pred(i,k) = 1.0;
      end
    end
   case 'ddATP'
    for i = 1:length( sequence )-1
      if ( sequence(i+1) == 'T' | sequence(i+1) == 'U' ); area_pred(i,k) = 1.0; end;
    end
   case 'ddCTP'
    for i = 1:length( sequence )-1
      if ( sequence(i+1) == 'G' ); area_pred(i,k) = 1.0; end;
    end
   case 'ddGTP'
    for i = 1:length( sequence )-1
      if ( sequence(i+1) == 'C' ); area_pred(i,k) = 1.0; end;
    end
   case 'ddTTP'
    for i = 1:length( sequence )-1
      if ( sequence(i+1) == 'A' ); area_pred(i,k) = 1.0; end;
    end
  end

  % following is useful for mutate-and-map data sets.
  mutpos = [];
  if isnumeric( data_types{k} ) % may be a number reflecting the mutation position.
    mutpos = round( data_types{k} );
  else
    mutpos = str2num( data_types{k} );
  end
  if ~isempty( mutpos ) & ~isnan( mutpos )
    area_pred(:,k) = 0.0;
    seqidx = mutpos - offset;
    if ( seqidx <= length( sequence) & seqidx >= 1 )
      area_pred(seqidx ,k) = 1.0;
    end
    if ~exist( 'bps' ) bps = convert_structure_to_bps( structure ); end;
    partner_idx = find_partner( seqidx, bps );
    if ( partner_idx >=1 & partner_idx <= length( sequence ) ) area_pred( partner_idx,k ) = 1.0; end;
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function partner_idx = find_partner( seqidx, bps )

partner_idx = 0;

for k = 1:size( bps, 1 );

  if bps(k,1) == seqidx; 
    partner_idx = bps(k,2); return;
  end

  if bps(k,2) == seqidx; 
    partner_idx = bps(k,1); return;
  end

end