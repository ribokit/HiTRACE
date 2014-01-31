function new_rdat =  cat_rdat_files( new_file, rdat_files );
%
%   cat_rdat_files( new_file, rdat_files );
%

if nargin==0; help( mfilename ); return; end;

for i = 1:length( rdat_files )
  rdats{i} = read_rdat_file( rdat_files{i} );
end
new_rdat = rdats{1};

% It is currently necessary that a lot of the basic information 
% between files is shared. Throw an error otherwise
for i = 2:length( rdats)
  if ~strcmp( rdats{i}.version, rdats{1}.version )
    fprintf( 'mismatch in version of rdat 1 and %d\n', i ); return;
  elseif ~strcmp( rdats{i}.name, rdats{1}.name )
      fprintf( 'mismatch in name of rdat 1 and %d\n', i ); return;
  elseif ~strcmp( rdats{i}.sequence, rdats{1}.sequence )
      fprintf( 'mismatch in sequence of rdat 1 and %d\n', i ); return;
  elseif ~strcmp( rdats{i}.structure, rdats{1}.structure )
      fprintf( 'mismatch in structure of rdat 1 and %d\n', i ); return;
  elseif length( rdats{i}.seqpos ) ~=  length( rdats{1}.seqpos ) 
      fprintf( 'mismatch in seqpos of rdat 1 and %d\n', i ); return;
  elseif rdats{i}.offset ~=  rdats{1}.offset
      fprintf( 'mismatch in offset of rdat 1 and %d\n', i ); return;
  end
end

% find all shared annotations.
new_annotations = rdats{1}.annotations;

for i = 2:length( rdats )
  rdat = rdats{i};

  match = zeros( 1, length( new_annotations ) );
  for j = 1:length( new_annotations )
    for k = 1:length( rdat.annotations )
      if strcmp( rdat.annotations{k}, new_annotations{j} )
	match( j ) = 1;
      end
    end
  end
  new_annotations = new_annotations( find( match ) );
  
end

%what annotations are specific to files?

for i = 1:length( rdats )
  annotations_to_add{ i } = {};

  rdat = rdats{i};
  match = zeros( 1, length( rdat.annotations ) );
  for k = 1:length( rdat.annotations )
    for j = 1:length( new_annotations )
      if strcmp( rdat.annotations{k}, new_annotations{j} )
	match( k ) = 1;
      end
    end
  end

  annotations_to_add{i} = rdat.annotations( find( ~match ) ) ;
  
end


new_rdat = rdats{1};
new_rdat.annotations = new_annotations;


if length( new_rdat.xsel ) > 0 | length( new_rdat.xsel_refine) > 0 | length( new_rdat.trace ) > 0 
  fprintf( 'WARNING! cat_rdat_files will not concatenate xsel, xsel_refine, or trace!\n');
end
new_rdat.xsel = [];
new_rdat.xsel_refine = [];
new_rdat.trace = [];


for i = 2:length( rdats )
  new_rdat.mutpos = [ new_rdat.mutpos, rdats{i}.mutpos ];
  new_rdat.comments = [ new_rdat.comments, rdats{i}.comments ];
  new_rdat.reactivity = [new_rdat.reactivity, rdats{i}.reactivity ];
  new_rdat.reactivity_error = [new_rdat.reactivity_error, rdats{i}.reactivity_error ];
end

count = 0;
for i = 1:length( rdats )
  rdat = rdats{i};
  for j = 1:length( rdat.data_annotations )
    count = count + 1;
    new_rdat.data_annotations{count} = [ rdat.data_annotations{j}, annotations_to_add{i} ];
  end
end

output_rdat_to_file( new_file,  new_rdat );

show_rdat( new_rdat );

