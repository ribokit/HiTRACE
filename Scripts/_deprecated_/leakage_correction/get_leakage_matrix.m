function lm = get_leakage_matrix( dye_names_full );

% where are we?
dirname = which( 'get_leakage_matrix.m' );
dirname = strrep( dirname, 'get_leakage_matrix.m', '' );

N = length( dye_names_full );
for i = 1:N

  lm(i,:) = zeros( 1, N );
  lm(i,i) = 1.0;

  if length( dye_names_full{i} ) == 0; continue;  end

  correction_file = [dirname, 'database/', dye_names_full{i}, '.txt' ];

  if ~exist( correction_file, 'file' )
    fprintf( 'Did not recognize %s! You must choose from the following:\n', dye_names_full{i} );
    dir( [dirname, 'database/*.txt'] );
    fprintf( ' ... or stick a new file in %s.\n', dirname );    
  end

  correction = textread( correction_file, '%f' );
  lm(i,:) = correction';
end


