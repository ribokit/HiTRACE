function area_pred = get_area_pred( sequence, data_types, structure );
% GET_PRED
%
% [ marks, area_pred, mutpos] = ( sequence, data_types, structure );
%
% Returns information on where bands should show up, based on input 
%  sequence, structure, and what modification/ladder reaction is in each lane.
%
% sequence  = full sequence
% data_types = cell of tags of modification reactions in each trace, e.g., 
%               {'SHAPE','SHAPE','ddTTP'}
% structure = structure in dot/bracket notation [give as '' if unknown]
%
% OUTPUT:
% area_pred = matrix of zeros and ones corresponding to where bands will 
%              show up.
%
% (C) R. Das, 2010-2013.
%

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
end


