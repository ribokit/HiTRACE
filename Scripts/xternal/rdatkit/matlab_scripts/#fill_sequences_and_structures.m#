function rdat = fill_sequences_and_structures( rdat );

% let's try to fill the "sequences" field if it isn't there.
rdat = fill_sequences_if_empty( rdat );
rdat = fill_structures_if_empty( rdat );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rdat = fill_sequences_if_empty( rdat );

if length( rdat.data_annotations ) == 0; return; end;

for i = 1:size(  rdat.reactivity, 2 )
  if ( i > length( rdat.sequences )  |  length(rdat.sequences{i} ) == 0 ) & i <= length( rdat.data_annotations )
    rdat.sequences{i} = rdat.sequence;
    
    data_annotation = rdat.data_annotations{i};
    for m = 1:length( data_annotation )

      c = str2cell( data_annotation{m},':' );

      if length(c) > 0  & strcmp( c{1}, 'sequence' )
	rdat.sequences{i} = c{2};
	continue;
      end
      
      if length(c) > 0  & strcmp( c{1}, 'mutation' )
      
	start_seq = '';
	mut_seq = '';
	mut_num = '';
	tag = cell2str( c(2:end), ':' );
      
	q = 1;	
	while ( q <= length(tag)  & isempty( str2num( tag(q) ) ) & tag(q)~='(' )
	  start_seq = [start_seq, tag(q) ];
	  q = q+1;
	end
	if q <= length(tag) & tag(q) == '('; q = q+1; end;
	while ( q <= length( tag ) &  (~isempty( str2num( tag(q) ) ) | tag(q)==':')  )
	  mut_num = [mut_num, tag(q) ];
	  q = q+1;
	end
	if q <= length(tag) & tag(q) == ')'; q = q+1; end;
	while q <= length( tag );
	  mut_seq = [mut_seq, tag(q) ];
	  q = q+1;
	end
	
	if  length( mut_num ) == 0 
	  if ~strcmp( start_seq, 'WT' )
	    fprintf( 'WARNING! Could not find mutation position in mutation annotation: %s\n', tag );
	  end
	  continue;
	end
	if  length( start_seq ) == 0 
	  fprintf( 'WARNING! Could not find starting nucleotide in mutation annotation: %s\n', tag );
	end
	if  length( mut_seq ) == 0 
	  fprintf( 'WARNING! Could not find mutation nucleotide in mutation annotation: %s\n', tag );
	end
	
	mutpos = str2num( mut_num ) - rdat.offset;
	if ( ~strcmp( rdat.sequence( mutpos ), start_seq ) )
	  mutpos = str2num( mut_num ); % perhaps specified without offset...
	    fprintf( 'WARNING! Mismatch between mutation nucleotides: %s (mutation annotation) vs. %s. (sequence)\n', rdat.sequence(mutpos), start_seq );
	  if ( strcmp( rdat.sequence( mutpos ), start_seq ) )
	    fprintf( 'OK, specified mutpos without taking into account offset...\n' );
	  end
	end
	rdat.sequences{i} = [rdat.sequence(1: (min(mutpos)-1) ), mut_seq, rdat.sequence( max(mutpos)+1 : end )];
      end
      
    end
  end
end
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rdat = fill_structures_if_empty( rdat );

if size( rdat.reactivity, 2 ) == 0; return; end;

for i = 1:size(  rdat.reactivity, 2 )

  if i > length( rdat.structures )  |  length(rdat.structures{i} ) == 0
    rdat.structures{i} = rdat.structure;
  
    data_annotation = rdat.data_annotations{i};
    for m = 1:length( data_annotation )

      c = str2cell( data_annotation{m},':' );
      if length(c) > 0  & strcmp( c{1}, 'structure' )
	rdat.structures{i} = c{2};
	continue;
      end
    end
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function c = str2cell(s, delim)

if ~exist( 'delim' ) 
  tabchar = sprintf( '\t' );
  if ~isempty(strfind( s, tabchar ) ) 
    delim = tabchar;
  else
    delim = ' '; 
  end
end;

rest = s;
i = 1;
c = {};
while length(rest)
    [t, rest] = strtok(rest, delim);
    c{i} = t;
    i = i + 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = cell2str( c, delim )
if  ~exist('delim','var'); delim = ' '; end;
s = '';
if length(c) > 0
s = c{1};
for i = 2:length(c); s = [s,delim,c{i}]; end
end
