function [ marks, area_pred, mutpos] = get_predicted_marks_SHAPE_DMS_CMCT( structure, sequence, offset,seqpos,data_type);
% GET_PREDICTED_MARKS_SHAPE_DMS_CMCT
%
% [ marks, area_pred, mutpos] = get_predicted_marks_SHAPE_DMS_CMCT( structure, sequence, offset,seqpos,data_type);
%
% Returns information on where bands should show up, based on input 
%  sequence, structure, and what modification/ladder reaction is in each lane.
%
% structure = structure in dot/bracket notation [give as '' if unknown]
% sequence  = full sequence
% offset    = integer to add to sequence position to get 'conventional' numbering
% seqpos    = numbering of nucleotides which give signals in capillary trace
% data_type = cell of tags of modification reactions in each trace, e.g., 
%               {'SHAPE','SHAPE','ddTTP'}
%
%
% Output: 
% marks = pairs of (mutpos,seqpos) where bands should show up. This format 
%         dates to the original use of hitrace to annotate mutate-and-map data, 
%         where mutpos is the sequence position that was mutated. More
%         generally mutpos contains a dummy integer that encodes the
%         kind of reaction/ladder, as follows
%             99999990. nomod
%             99999991. SHAPE 
%             99999992. DMS
%             99999993. CMCT
%             99999994. A (ddTTP ladder)
%             99999995. C (ddGTP ladder)
%             99999996. G (ddCTP ladder)
%             99999997. U (ddATP ladder)
%             99999998. UV (double-pyrimidine ladder)
%         This may be deprecated in later versions of hitrace.
% area_pred = matrix of zeros and ones corresponding to where bands will 
%              show up.
% mutpos = the integer 'mutpos' code for each capillary trace. This may be
%           deprecated in later versions of hitrace.
%
%
% (C) R. Das, 2010-2013.
%
FMN_SHAPE = [0 0 0 1 0 0 0 0 0 0 0];
FMN_DMS = [0 0 0 1 0 1 0 1 1 0 0];

marks = [];
area_pred = [];
mutpos = [];

MUTPOS_OFFSET = 99999990;

if isempty('seqpos')
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

% sequencing ladders
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


% UV -- pyrimidines.
for i = 2:length( structure )
  if ( sequence(i) == 'C'  | sequence(i) == 'U'  ) & ...
	( sequence(i-1) == 'C'  | sequence(i-1) == 'U'  ) 
    marks = [ marks; MUTPOS_OFFSET+8, i+offset ];
  end
end

d_pred = zeros( length(sequence), 8 );
for i = [1:8]
  d_pred( marks( find( marks(:,1) == MUTPOS_OFFSET+i ), 2 ) - offset, i ) = 1.0;
end

% 'ideal' predicted data -- needed for background subtraction.
% probably should tuck into a function...

area_pred = [];
for k = 1:length( mutpos )
  area_pred(:,k) = zeros( length(seqpos),1 );
  if mutpos(k) > MUTPOS_OFFSET    
    area_pred(:,k) = d_pred( seqpos-offset, mutpos(k)-MUTPOS_OFFSET );
    if mutpos(k) == MUTPOS_OFFSET + 1
        reaction_target = find(structure == 'F');
        if length(reaction_target) == length(FMN_SHAPE)
            area_pred(seqpos(reaction_target),k) = FMN_SHAPE';
        end
    elseif mutpos(k) == MUTPOS_OFFSET + 2
        reaction_target = find(structure == 'F');
        if length(reaction_target) == length(FMN_DMS)
            area_pred(seqpos(reaction_target),k) = FMN_DMS';
        end
    end
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

OK_labels = {'SHAPE','DMS','CMCT','ddTTP','ddGTP','ddCTP','ddATP','UV','nomod'};

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