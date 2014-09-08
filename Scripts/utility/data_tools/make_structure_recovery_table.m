function [ num_native_stems, all_num_native_stems_recovered, all_num_false_prediction_stems ] = make_structure_recovery_table( tags, all_native_structure, all_structures, rna_structure_run_type );
	   
%function [native_stems, all_stems, all_correspondence_to_native ] = make_structure_recovery_table( tags, all_native_structure, all_structures, rna_structure_run_type );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analyze accuracies.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:length( tags )
  native_structure = all_native_structure{i};
  native_stems = parse_stems( native_structure );

  %tail_pos = figure_out_tails( native_structure );
  [crystpos, tail_pos ] = get_tail_crystpos( tags{i}, native_structure );

  num_native_stems = length( native_stems );
  all_num_native_stems(i) = num_native_stems;

  all_num_extra_bps_added_to_stem{i} = [];
  all_num_missing_bps_in_stem{i}     = [];
  all_stems{i}                       = {};
  all_correspondence_to_native{i}    = {};

  structures       = all_structures{i};
  tags{i}
  structures

  inpath_common = 'OUTPUT1D/'; 
  seq_file = [inpath_common,tags{i},'_sequence.seq'];
  sequence = read_seq_file( seq_file );

  num_res_no_tails(i) = length(all_native_structure{i})-length(tail_pos);
  
  for j = 1:length( structures )
    structure = structures{j};
    stems = parse_stems( structure );
    stems = filter_stems( stems, tail_pos );

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Native recovery
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [ stem_correspondence, num_extra_bps_added_to_stem, num_missing_bps_in_stem ] = find_shared_stems( native_stems, stems );
      
    all_num_extra_bps_added_to_stem{i}{j} =  num_extra_bps_added_to_stem;
    all_num_missing_bps_in_stem{i}{j}     =  num_missing_bps_in_stem;

    num_native_stems_recovered = length( setdiff( unique( stem_correspondence ), [0] ) );    
    num_false_prediction_stems = length( stems ) - num_native_stems_recovered;
      
    all_num_native_stems_recovered(i,j) = num_native_stems_recovered;
    all_num_false_prediction_stems(i,j) = num_false_prediction_stems;

    [ num_native_bps(i,j), num_recovered_bps(i,j), num_extra_bps(i,j), num_extra_AU_bps(i,j), num_extra_GC_bps(i,j), num_missed_bps(i,j), num_missed_AU_bps(i,j), num_missed_GC_bps(i,j), num_false_bps(i,j), num_recovered_bps_allow_shift(i,j), num_false_bps_allow_shift(i,j)] = get_fine_error_info( stem_correspondence, native_stems, stems, sequence );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties of stems
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    all_stems{i}{j} = stems;
    
    correspondence_to_native = find_shared_stems( stems, native_stems );
    all_correspondence_to_native{i}{j} = correspondence_to_native;

  end
  
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% within recovered helices, what errors are there?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


table_file = 'bp_breakdown_table.txt' ;
fid_table = fopen( table_file, 'w' );

for fid = [1 fid_table ]

  j = length( structures );
  num_extra_GU_bps =  num_extra_bps - num_extra_AU_bps - num_extra_GC_bps;
  num_missed_GU_bps =  num_missed_bps - num_missed_AU_bps - num_missed_GC_bps;
  fprintf( fid, '%20s                 missed          extra   \n','' );
  fprintf( fid, '%20s   tot   rec   A-U G-U G-C   A-U G-U G-C\n','RNA' );
  for i = 1:length( tags )
    fprintf( fid, '%20s   %3d   %3d   %3d %3d %3d   %3d %3d %3d\n', tags{i}, ...
	     num_native_bps(i,j), num_recovered_bps(i,j), ...
	     num_missed_AU_bps(i,j), num_missed_GU_bps(i,j), num_missed_GC_bps(i,j), ...
	     num_extra_AU_bps(i,j),  num_extra_GU_bps(i,j), num_extra_GC_bps(i,j)  );
  end
  fprintf( fid, '%20s   %3d   %3d   %3d %3d %3d   %3d %3d %3d\n', 'Total', ...
	   sum(num_native_bps(:,j)), sum(num_recovered_bps(:,j)), ...
	   sum(num_missed_AU_bps(:,j)), sum(num_missed_GU_bps(:,j)), sum(num_missed_GC_bps(:,j)),...  
	   sum(num_extra_AU_bps(:,j)),  sum(num_extra_GU_bps(:,j)), sum(num_extra_GC_bps(:,j)) );
  fprintf( '\n');
end

fclose( fid_table );
fprintf( '\nCreated table file: %s\n\n', table_file );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make (fine) table of errors in 'correct' stems
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bins = [0:3];
count_num_extra_bps_added_to_stem = zeros( length(tags), length(bins) );
count_num_missing_bps_in_stem = zeros( length(tags), length(bins) );
for i = 1:length( tags )
  j = length( structures );

  num_extra_bps_added_to_stem = all_num_extra_bps_added_to_stem{i}{j};
  for k = 1:length( num_extra_bps_added_to_stem )
    bin = find( bins == num_extra_bps_added_to_stem( k ) );
    if ~isempty( bin ); count_num_extra_bps_added_to_stem( i, bin ) = count_num_extra_bps_added_to_stem( i, bin ) + 1; end;
  end

  num_missing_bps_in_stem = all_num_missing_bps_in_stem{i}{j};
  for k = 1:length( num_missing_bps_in_stem )
    bin = find( bins == num_missing_bps_in_stem( k ) );
    if ~isempty( bin ); count_num_missing_bps_in_stem( i, bin ) = count_num_missing_bps_in_stem( i, bin ) + 1; end;
  end
end

fprintf( '\n' );
fprintf( '%20s   num added     num missing\n', '' );
fprintf( '%20s', 'RNA' );
for bin = 1:length( bins )
  fprintf( ' +%d', bins(bin) )
end
  fprintf( '    ' );
for bin = 1:length( bins )
  fprintf( ' -%d', bins(bin) )
end
fprintf( '\n' );

for i = 1:length( tags )
  fprintf( '%20s', tags{i} );
  for bin = 1:length( bins )
    fprintf( ' %2d', count_num_extra_bps_added_to_stem(i,bin ) );
  end
  fprintf( '    ' );
  for bin = 1:length( bins )
    fprintf( ' %2d', count_num_missing_bps_in_stem(i,bin ) );
  end
  fprintf( '\n' );
end
fprintf( '%20s', 'Total' );
for bin = 1:length( bins )
  fprintf( ' %2d', sum(count_num_extra_bps_added_to_stem(:,bin )) );
end
fprintf( '    ' );
for bin = 1:length( bins )
  fprintf( ' %2d', sum(count_num_missing_bps_in_stem(:,bin )) );
end
fprintf( '\n' );

fprintf( '\n' );



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make bp-level table... allow +/-1 shift
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
table_file = 'bp_recovery_table_allow_shift.txt' ;
print_table( table_file, structures, rna_structure_run_type, tags, num_res_no_tails, num_native_bps(:,1), num_recovered_bps_allow_shift, num_false_bps_allow_shift );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make bp-level table.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
table_file = 'bp_recovery_table.txt' ;
print_table( table_file, structures, rna_structure_run_type, tags, num_res_no_tails, num_native_bps(:,1), num_recovered_bps, num_false_bps );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make table.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
table_file = 'structure_recovery_table.txt' ;
print_table( table_file, structures, rna_structure_run_type, tags, num_res_no_tails, all_num_native_stems, all_num_native_stems_recovered, all_num_false_prediction_stems );



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function    [ num_native_bps, num_recovered_bps, num_extra_bps, num_extra_AU_bps, num_extra_GC_bps, num_missed_bps, num_missed_AU_bps, num_missed_GC_bps, num_false_bps, num_recovered_bps_allow_shift, num_false_bps_allow_shift ] = get_fine_error_info( stem_correspondence, native_stems, stems, sequence );

num_native_bps = 0;

num_recovered_bps = 0;
num_recovered_bps_allow_shift = 0;

num_extra_bps = 0;
num_extra_AU_bps = 0;
num_extra_GC_bps = 0;

num_missed_bps = 0;
num_missed_AU_bps = 0;
num_missed_GC_bps = 0;

for i = 1:length( native_stems )

  native_stem = native_stems{i};
  num_native_bps = num_native_bps + size( native_stem, 1 );

  j = stem_correspondence( i );

  if  j > 0
    stem = stems{j};
    
    ok = 0;
    ok_allow_shift = 0;
    
    found_match_in_stem = zeros( 1, size( stem, 1 ) );
    
    for p = 1:size( native_stem, 1 )

      for q = 1:size( stem, 1)
	
	if ( abs( native_stem( p, 1) - stem( q, 1 ) ) <= 1 & ...
	     abs( native_stem( p, 2) - stem( q, 2 ) ) <= 1 )
	  ok_allow_shift = 1;
	end
	
	if ( native_stem( p, 1) == stem( q, 1 ) & ...
	     native_stem( p, 2) == stem( q, 2 ) )
	  ok = 1;
	  found_match_in_stem(q) = 1; break;
	end
      end
          
      if  ok
	num_recovered_bps = num_recovered_bps + 1;
      else
	num_missed_bps = num_missed_bps + 1;

	r = native_stem(p,1); s = native_stem(p,2);
	if ( sequence(r) == 'A' &  sequence(s) == 'U' ) | ...
	      ( sequence(s) == 'A' &  sequence(r) == 'U' ) 
	  num_missed_AU_bps = num_missed_AU_bps + 1;
	elseif ( sequence(r) == 'G' &  sequence(s) == 'C' ) | ...
	      ( sequence(s) == 'C' &  sequence(r) == 'G' ) 
	  num_missed_GC_bps = num_missed_GC_bps + 1;	
	end
      end
      
    end

    for q = 1:size( stem, 1)
      if ~found_match_in_stem( q) 
	num_extra_bps = num_extra_bps + 1;
	r = stem(q,1); s = stem(q,2);
	if ( sequence(r) == 'A' &  sequence(s) == 'U' ) | ...
	      ( sequence(s) == 'A' &  sequence(r) == 'U' ) 
	  num_extra_AU_bps = num_extra_AU_bps + 1;
	elseif ( sequence(r) == 'G' &  sequence(s) == 'C' ) | ...
	      ( sequence(s) == 'C' &  sequence(r) == 'G' ) 
	  num_extra_GC_bps = num_extra_GC_bps + 1;	
	end	
      end
    end
  end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_recovered_bps_allow_shift = 0;

for i = 1:length( native_stems )

  native_stem = native_stems{i};

  for q = 1:size( native_stem, 1)

    ok_allow_shift = 0;

    for r = 1:length( stems )
      stem = stems{r};
      
      for p = 1:size( stem, 1 )
	
	if ( abs( stem( p, 1) - native_stem( q, 1 ) ) <= 1 & ...
	     abs( stem( p, 2) - native_stem( q, 2 ) ) <= 1 )
	  ok_allow_shift = 1; break;
	end
	
      end
      
      if ( ok_allow_shift == 1); break; end;
      
    end

    if  ok_allow_shift
      num_recovered_bps_allow_shift = num_recovered_bps_allow_shift + 1;
    end

  end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_false_bps = 0;
num_false_bps_allow_shift = 0;

for i = 1:length( stems )

  stem = stems{i};

  for q = 1:size( stem, 1)

    ok = 0;
    ok_allow_shift = 0;

    for r = 1:length( native_stems )
      native_stem = native_stems{r};
    
      for p = 1:size( native_stem, 1 )

	if ( abs( native_stem( p, 1) - stem( q, 1 ) ) <= 1 & ...
	     abs( native_stem( p, 2) - stem( q, 2 ) ) <= 1 )
	  ok_allow_shift = 1;
	end

	if ( native_stem( p, 1) == stem( q, 1 ) & ...
	     native_stem( p, 2) == stem( q, 2 ) )
	  ok = 1;
	end
      end
      
      if ( ok == 1); break; end;
      
    end

    if  ~ok
      num_false_bps = num_false_bps + 1;
    end
    if  ~ok_allow_shift
      num_false_bps_allow_shift = num_false_bps_allow_shift + 1;
    end

  end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function print_table( table_file, structures, rna_structure_run_type, tags, num_res_no_tails, all_num_native_bps, all_num_native_bps_recovered, all_num_false_prediction_bps );

fid_table = fopen( table_file, 'w' );

for fid = [1 fid_table ]

  fprintf( fid, '\n' );
  fprintf( fid,   '%20s %4s  %4s ',' ',' ',' ' );
  for j = 1:length( structures )
    fprintf( fid, '  %10s', rna_structure_run_type{j} );
  end
  fprintf( fid, '\n');
  
  fprintf( fid,   '%20s %4s  %4s', 'RNA','LEN','TOT' );
  for j = 1:length( structures )
    fprintf( fid, '     %3s  %3s', 'TP','FP' );
  end
  fprintf( fid, '\n');
  
  for i = 1:length( tags )
    fprintf( fid, '%20s %4d  %4d', tags{i}, num_res_no_tails(i), all_num_native_bps(i) );
    for j = 1:length( structures )
    fprintf( fid, '     %3d  %3d', all_num_native_bps_recovered(i,j), all_num_false_prediction_bps(i,j) );
    end
    fprintf( fid, '\n' );
  end

  fprintf( fid, '%20s  %3d   %3d', 'Total', sum( num_res_no_tails ), sum( all_num_native_bps(:,1) ) );
  for j = 1:length( structures )
    fprintf( fid, '     %3d  %3d', sum(all_num_native_bps_recovered(:,j)), sum(all_num_false_prediction_bps(:,j)) );
  end
  fprintf( fid, '\n' );
  fprintf( fid, '-------------------------------------------------------------\n');
  fprintf( fid, '%20s %4s   %2s','FNR',' ',' ');
  for j = 1:length( structures )
    TP(j) = sum( all_num_native_bps_recovered(:,j) );
    FN(j) = sum( all_num_native_bps ) - TP(j);
  fprintf( fid, '   %10.1f', 100*FN(j)/(TP(j)+FN(j)) );
  end
  fprintf( fid, '\n' );
  fprintf( fid, '%20s %4s   %2s','FDR',' ',' ');
  for j = 1:length( structures )
    FP(j) = sum( all_num_false_prediction_bps(:,j) );
    fprintf( fid, '   %10.1f%',100* FP(j)/(TP(j)+FP(j)));
  end
  fprintf( fid, '\n' );  
  fprintf( fid, '-------------------------------------------------------------\n');

  fprintf( fid, '%20s %4s   %2s','Sensitivity',' ',' ');
  for j = 1:length( structures )
    fprintf( fid, '  %10.1f', 100*TP(j)/(TP(j)+FN(j)) );
  end
  fprintf( fid, '\n' );
  fprintf( fid, '%20s %4s   %2s','PPV',' ',' ');
  for j = 1:length( structures )
    fprintf( fid, '  %10.1f%',100* TP(j)/(TP(j)+FP(j)));
  end
  fprintf( fid, '\n' );  

end

fclose( fid_table );
fprintf( '\nCreated table file: %s\n\n', table_file );
