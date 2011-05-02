function [ marks, area_pred, mutpos] = get_predicted_marks_SHAPE_DMS_CMCT( structure, sequence, offset,seqpos,data_type);
% 99999990. nomod
% 99999991. SHAPE 
% 99999992. DMS
% 99999993. CMCT
% 99999994. A (ddTTP ladder)
% 99999995. C (ddGTP ladder)
% 99999996. G (ddCTP ladder)
% 99999997. U (ddATP ladder)

marks = [];
area_pred = [];
mutpos = [];

MUTPOS_OFFSET = 99999990;

if ~isempty('seqpos')
    seqpos = length(sequence)-20 - [1:(length(sequence)-20)] + 1 + offset;
end

if iscell( data_type )
  mutpos = figure_out_mutpos( data_type, MUTPOS_OFFSET );
  
  ok_pos = find( ~isnan( mutpos ) & mutpos < MUTPOS_OFFSET );
  for k = unique(mutpos(ok_pos));  marks = [ marks;  k,k ]; end; 
else
  mutpos = data_type;
end


% SHAPE
if length( structure ) == 0; 
  structure = sequence; 
  for i = 1:length( structure ); structure(i) = '.'; end;
end

for i = 1:length( structure )
  if structure(i) == '.'
    marks = [ marks; MUTPOS_OFFSET+1, i+offset ];
  end
end


% DMS
for i = 1:length( structure )
  if structure(i) == '.' & ( sequence(i) == 'A'  | sequence(i) == 'C'  ) 
    marks = [ marks; MUTPOS_OFFSET+2, i+offset ];
  end
end

% CMCT
for i = 1:length( structure )
  if structure(i) == '.' & ( sequence(i) == 'G'  | sequence(i) == 'U'  | sequence(i) == 'T'  ) 
    marks = [ marks; MUTPOS_OFFSET+3, i+offset ];
  end
end

% CMCT
for i = 1:length( sequence )-1

  switch sequence(i+1)
   case 'A'
    marks = [ marks; MUTPOS_OFFSET+4, i+offset ];
   case 'C'
    marks = [ marks; MUTPOS_OFFSET+5, i+offset ];
   case 'G'
    marks = [ marks; MUTPOS_OFFSET+6, i+offset ];
   case {'T','U'}
    marks = [ marks; MUTPOS_OFFSET+7, i+offset ];
  end
  
end


d_pred = zeros( length(sequence), 7 );
for i = [1:7]
  d_pred( marks( find( marks(:,1) == MUTPOS_OFFSET+i ), 2 ) - offset, i ) = 1.0;
end


% 'ideal' predicted data -- needed for background subtraction.
% probably should tuck into a function...

area_pred = [];
for k = 1:length( mutpos )
  area_pred(:,k) = zeros( length(seqpos),1 );
  if mutpos(k) > MUTPOS_OFFSET    
    area_pred(:,k) = d_pred( seqpos-offset, mutpos(k)-MUTPOS_OFFSET );
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  mutpos = figure_out_mutpos( data_type, MUTPOS_OFFSET );
% 1. SHAPE 
% 2. DMS
% 3. CMCT
% 4. A (ddTTP ladder)
% 5. C (ddGTP ladder)
% 6. G (ddCTP ladder)
% 7. U (ddATP ladder)

mutpos = [];

OK_labels = {'SHAPE','DMS','CMCT','ddTTP','ddGTP','ddCTP','ddATP','nomod'};

for i = 1:length( data_type )

  mutpos(i) = NaN;

  for j = 1: length( OK_labels )
    if ~isempty( strfind( data_type{i}, OK_labels{j} )  )
      mutpos( i ) = MUTPOS_OFFSET + j;
      break;
    end
  end

  
  if isnan( mutpos(i) )
    possible_number = str2num( data_type{i} );
    if ~isempty( possible_number )
      mutpos(i) = possible_number;
    else
      fprintf( 'Unrecognized label type, %s\n', data_type{i} );
      fprintf( 'OK labels are: ')
      char(OK_labels)
      return;
    end
  end

  if ( mutpos(i) == MUTPOS_OFFSET + length( OK_labels ) );    mutpos(i) = MUTPOS_OFFSET; end;
  
end