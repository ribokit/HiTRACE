function output_rdat_to_file( filename, rdat );
%
% output_rdat_to_file( filename, rdat );
%
% Copyright R. Das, P. Cordero, Stanford University, 2010,2011
%

if nargin == 0; help( mfilename ); return; end;

fprintf( 'About to create file: %s\n', filename );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare text for output.

s = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%s = [s, 'RDAT_VERSION 0.22\n'];
%s = [s, 'RDAT_VERSION 0.23\n']; % use sorted sequence order; spacer in SEQPOS
%s = [s, 'RDAT_VERSION 0.24\n']; % use REACTIVITY instead of AREA_PEAK
%s = [s, 'RDAT_VERSION\t0.30\n']; % use tab-delimiters!
%s = [s, 'RDAT_VERSION\t0.31\n']; % reorder so that SEQPOS & OFFSET occur closer to data
%s = [s, 'RDAT_VERSION\t0.32\n']; % provide nucleotides in SEQPOS
s = [s, 'RDAT_VERSION\t0.33\n']; % get rid of MUTPOS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s = [s, 'NAME\t', rdat.name,'\n'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s = [s, 'SEQUENCE\t', rdat.sequence, '\n'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% its nice to put structure right after sequence, since it should be the same length.
if length(rdat.structure ) > 0
  s = [s, 'STRUCTURE\t', rdat.structure, '\n'];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length( rdat.comments ) > 0
  s = [s,'\n'];
  for i = 1:length( rdat.comments );  s = [s, 'COMMENT\t', rdat.comments{i},'\n']; end;
  s = [s,'\n'];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length( rdat.annotations )>0;    
  s = [s, 'ANNOTATION\t', cell2str( rdat.annotations,'\t'), '\n']; 
end

fid = fopen(filename, 'w');
fprintf( fid, s );
s = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this needs to happen line by line -- running out of memory actually.
if  ~isempty( rdat.data_annotations ) 

  print_data_annotations = 0;
  for i=1:length( rdat.data_annotations )
    if length( rdat.data_annotations{i} ) > 0;
      print_data_annotations = 1; break;
    end
  end  

  if print_data_annotations
    %s = [s,'\n'];
    fprintf( fid, '\n' );
    for i=1:length( rdat.data_annotations )
      %s = [s, 'ANNOTATION_DATA:', int2str_exact(i,6),'        ',cell2str( rdat.data_annotations{i},' '),'\n'];
      fprintf( fid, ['ANNOTATION_DATA:', int2str_exact(i,6),'\t',cell2str( rdat.data_annotations{i},'\t'),'\n'] );
    end
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf( fid, '\n' );
fprintf( fid,  ['OFFSET\t', num2str( rdat.offset, '\t%d'),'\n'] );
fprintf( fid, '\n' );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf( fid,  ['SEQPOS\t']);
for i = 1:length( rdat.seqpos )
  m = rdat.seqpos(i) - rdat.offset;
  if ( m < 1 | m > length( rdat.sequence ) )
    fprintf( 'Error in SEQPOS %d -- please check OFFSET and SEQUENCE\n', rdat.seqpos(i) );
    error('Error in seqpos');
  end
  fprintf( fid, '\t%s%d', rdat.sequence(m), rdat.seqpos(i) );
end
fprintf( fid, '\n' );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%if length( rdat.mutpos ) > 0
%  s = [s, 'MUTPOS\t', strrep( num2str(rdat.mutpos,'\t%d'), 'NaN', ' WT' ), '\n'];
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OK, the good stuff.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%s = [s,'\n'];
%fprintf( fid, '\n' );
num_lanes = size( rdat.reactivity, 2  );
for i=1:num_lanes
  %s = [s, 'REACTIVITY:', int2str_exact(i,6), '              ', num2str( rdat.reactivity(:,i)', ' %9.4f'),'\n'];
   fprintf( fid, ['REACTIVITY:', int2str_exact(i,6), '\t', num2str( rdat.reactivity(:,i)', '\t%9.4f'),'\n' ]);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Should be specified for standard states.
if ~isempty( rdat.reactivity_error )
  %s = [s,'\n'];
  fprintf( fid, '\n' );
  for i=1:size( rdat.reactivity_error, 2 )
    %s = [s, 'REACTIVITY_ERROR:', int2str_exact(i,6), '        ', num2str( rdat.reactivity_error(:,i)', ' %9.4f'),'\n'];
    fprintf( fid, [ 'REACTIVITY_ERROR:', int2str_exact(i,6), '\t', num2str( rdat.reactivity_error(:,i)', '\t%9.4f'),'\n'] );
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'Raw data' -- optional
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty( rdat.xsel )
  s = [s,'\n'];
  if length( rdat.xsel ) ~= length( rdat.seqpos ); fprintf(  'mismatch in length between xsel and seqpos!\n' ); return; end;
  s = [s, 'XSEL\t', num2str( rdat.xsel, '\t%8.2f'), '\n'];
end

if ~isempty( rdat.xsel_refine )
  s = [s,'\n'];
  num_lanes = size( rdat.xsel_refine, 2 );
  for i=1:num_lanes
    s = [s,   'XSEL_REFINE:',  int2str_exact(i,6), '\t', num2str( rdat.xsel_refine(:,i)', '\t%8.2f'),'\n'];     
  end
  s = [s,'\n'];
end

if ~isempty( rdat.trace )
  Nprogress = 20;
  fprintf( 'Compiling trace: '); 
  progress = init_progress_bar( Nprogress ); 
  num_lanes = size( rdat.trace, 2 );

  s = [s,'\n'];
  for i=1:num_lanes
    s = [s, 'TRACE:',  int2str_exact(i,6), '\t', num2str( rdat.trace(:,i)', '\t%9.4f'),'\n'];    
    % update progress bar.
    progress = update_progress_bar( progress, Nprogress, i, num_lanes );
  end
  fprintf( '\n' );
end

%fprintf( 'Finished outputting text\n' );
% Print out the file.
fprintf(fid, s);
%fprintf( '\n' );
fclose( fid );

% quick check on names.
check_filename_vs_annotation( rdat, filename );

% quick check on any specified mutations
check_mutations_vs_sequence( rdat );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function ok = is_ok_data_type( data_type )
%ok_data_types = {'DMS','CMCT','SHAPE','nomod','ddATP','ddCTP','ddTTP','ddGTP'};
%ok = 0;
%for i = 1:length( ok_data_types )
%  if ( strcmp( data_type, ok_data_types{i} ) )
%    ok = 1; break
%  end
%end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = int2str_exact( i, n )
s = int2str( i );

for j = 1:(n-length(s))
  %s = [' ',s];
  s = [s,' '];  % for v0.23, where integer is attached to tag ('REACTIVITY') with a colon.
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = cell2str( c, delim )
s = '';
c = strrep( c, '%','%%'); % otherwise % is assumed to be an escape string.
if ~isempty(c)
  s = c{1};
  for i = 2:length(c); 
    s = [s,delim,c{i} ]; 
  end
end
s = char(s);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function progress = init_progress_bar( N )

fprintf( '[')
for j = 1:N; fprintf( ' ' );end;
fprintf(']' );
progress = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function progress = update_progress_bar( prev_progress, N, i, num_lanes )

progress = floor(N*i/num_lanes);
if ( progress ~= prev_progress )
  for j = 1:(N+1); fprintf( '\b' ); end; % remove previous progress bar.
  for j = 1:progress; fprintf( '=' );end;
  for j = (progress+1):N; fprintf( ' ');end;
  fprintf( ']' );
end    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function check_filename_vs_annotation( rdat, filename );

ok_data_types = {'DMS','CMCT','SHAPE','nomod','NOMOD'};

for j = 1:length( ok_data_types )
  if strfind( filename, ok_data_types{j} ) 
    for k = 1:length( ok_data_types )      
      if isempty(strfind( ok_data_types{k}, ok_data_types{j} ))
	for q = 1:length( rdat.annotations )
	  if ~isempty( strfind( rdat.annotations{q}, ok_data_types{k} ) ) 
	    fprintf( 'WARNING! filename has %s in name, but annotations includes modifier: %s\n', ok_data_types{j}, rdat.annotations{q} );
	  end
	end
      end    
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function check_mutations_vs_sequence( rdat );

for i = 1:length( rdat.data_annotations )
  data_annotation = rdat.data_annotations{i};

  for j = 1:length( data_annotation )
    mutation_annotation = data_annotation{j};
    [tag,mutation_string] = strtok( mutation_annotation, 'mutation:' );
    if length( mutation_string ) > 0 & ~strfind( mutation_string, 'WT' )
      [seqchar, num ] = parse_seqchar_number( mutation_string );
      if ( seqchar ~= rdat.sequence( num-rdat.offset ) ) 
	fprintf( 'WARNING! Mismatch in sequence! %c  vs.  %c', seqchar, mutation_string );
      end
    end
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [seqchar, num ] = parse_seqchar_number( construct_name );

% find the last number in the construct name.
in_number = 0;
num_start = 0;
num_end = 0;
num = 0;
seqchar = '?';

for i = length( construct_name ):-1:1
  if ( ~isempty( str2num( construct_name(i) ))  && ~in_number )
    num_end = i;
    in_number = 1;
  end
  if ( isempty( str2num( construct_name(i) ) ) && in_number )
    num_start = i+1; break;
  end
end      

if ( num_start<2 | num_end<3 ) 
  return;
end

num = str2num( construct_name( num_start:num_end ) );
seqchar = construct_name( num_start - 1 );