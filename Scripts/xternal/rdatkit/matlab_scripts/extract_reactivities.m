function [reactivities, seqpos ] = extract_reactivities( rdat,  outfile, idx );

if ischar( rdat ); rdat = read_rdat_file( rdat ); end;
if ~exist( 'idx','var') idx = 1; end;

fid = fopen( outfile, 'w' );
for i = 1:length( rdat.seqpos );
  fprintf( fid, '%d %9.5f\n', rdat.seqpos(i), rdat.reactivity(i,idx) );
end

fclose( fid );
fprintf( 'Created: %s\n', outfile );
