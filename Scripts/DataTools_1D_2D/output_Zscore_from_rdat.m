function [ Z, mutpos, seqplot ] = output_Zscore_from_rdat( outfile, rdat_files, d_nomod, MEAN_EXPOSURE_CUTOFF, ZSCORE_CUTOFF, APPLY_ZSCORE_CUTOFF, ONLY_A_C, print_stuff );

if ~exist( 'MEAN_EXPOSURE_CUTOFF' ) MEAN_EXPOSURE_CUTOFF = 1.0; end;
if ~exist( 'APPLY_ZSCORE_CUTOFF' ) APPLY_ZSCORE_CUTOFF = 1; end;
if ~exist( 'ZSCORE_CUTOFF' ) ZSCORE_CUTOFF = 0.0; end;
if ~exist( 'ONLY_A_C' ) ONLY_A_C = 0; end;
if ~exist( 'print_stuff' ); print_stuff = 0; end

Z = []; mutpos = [];seqplot = [];
if length( rdat_files ) == 0; return; end;

if ~iscell( rdat_files ); rdat_files = { rdat_files }; end;

NRES = 0;
if ~isempty( d_nomod); 
  if ischar( d_nomod )
    d_nomod = read_rdat_files( d_nomod);
  end
  NRES = length( d_nomod.sequence ); 
  Z_sum = zeros( NRES, NRES );
end;

mut_weights_sum = zeros( 1, NRES );

for i = 1:length( rdat_files )

  d = read_rdat_file( rdat_files{i} );

  if ( NRES == 0 )
    NRES = length(d.sequence);
    Z_sum = zeros( NRES, NRES );
    mut_weights_sum = zeros( 1, NRES );
  elseif ( length(d.sequence) ~= NRES ) 
    fprintf( 'WARNING! WARNING! NRES problem! %s\n', rdat_files{i} );
  end

  [ Z, mutpos ] = get_Zscore_and_apply_filter( d, d_nomod, MEAN_EXPOSURE_CUTOFF, ZSCORE_CUTOFF, APPLY_ZSCORE_CUTOFF, 0 );

  % just a check
  if ~isempty(strfind( rdat_files{i}, 'P4P6'));    Z( 176-d.offset, : ) = 0.0;   end;

  Z_sum = Z_sum + Z; % + smooth2d( Z );

  
  mut_weights = ( sum( Z, 1 ) < 0 );
  mut_weights_sum = mut_weights_sum + mut_weights;

end

for i = 1:NRES
  if ( mut_weights_sum(i) > 0 )
    Z(:,i) = Z_sum(:,i)/mut_weights_sum(i);    
  end
end

seqplot = [1:length( d.sequence )] + d.offset;
plot_and_save(Z, seqplot,  outfile, print_stuff );

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ Zscore_full, mutpos ] = get_Zscore_and_apply_filter( d, d_nomod, MEAN_EXPOSURE_CUTOFF, ZSCORE_CUTOFF, APPLY_ZSCORE_CUTOFF, ONLY_A_C  );

Zscore = [];

a_input = d.area_peak;
normbins = [ 10: size( a_input,1)-10 ];
a = quick_norm( a_input, normbins );

for i = 1:size( a, 1 )
  Zscore(i,:) =  ( a(i,:) - mean( a(i,:)) )/ std( a(i,:),0,2);
  mean_exposure(i) = mean( a(i,:) );
end

% special case -- what if some columns are different experiments.
diff_ex_cols = [];
for j = 1:size( a, 2 )
  data_annotation = d.data_annotations{j};
  for k = 1:length( data_annotation )
    if ( strfind( data_annotation{k}, 'differentExperiment' ) )
      diff_ex_cols = [ diff_ex_cols, j ];
      break;
    end
  end
end
%if ~isempty( 'diff_ex_cols' )
%  diff_ex_cols
%end

for i = 1:size( a, 1 )
  Zscore(i,diff_ex_cols) =  ( a(i,diff_ex_cols) - mean( a(i,diff_ex_cols)) )/ std( a(i,diff_ex_cols),0,2);
end

  
% Don't include rows that have uniformly high chemical modification --
% mutate/map looks for 'release' of protected residues.
for i = 1:size( a, 1 )
  if ( mean_exposure(i) > MEAN_EXPOSURE_CUTOFF  )
    Zscore(i,:) = 0.0;
  end
end

if ( ONLY_A_C )
  for i = 1:size( a, 1 )
    seqchar =   d.sequence( d.seqpos( i ) - d.offset );
    if ( seqchar ~='A' & seqchar ~='C' )
      Zscore(i,:) = 0.0;
    end  
  end
end

% Filter to make sure Z-score is above some cutoff.
if APPLY_ZSCORE_CUTOFF
  Zscore = max( ( Zscore - ZSCORE_CUTOFF ), 0.0 );
end

% filter for outlier columns (typically bad RNA in that well...)
%MEAN_Z_FILTER = 1;
%if MEAN_Z_FILTER
%  mean_Z = mean( abs( Zscore ) );
%
%  for j = 1:size( a_input, 2 );
%    if ( mean_Z(:,j) > mean(mean_Z) * 4 )
%      Zscore(:,j) = 0.0;
%    end
%  end
%end

for j = 1:size( a, 2 )
  data_annotation = d.data_annotations{j};
  for k = 1:length( data_annotation )
    if ( strfind( data_annotation{k}, 'badQuality' ) )
      Zscore(:,j) = 0.0;
    end
  end
end


% If 'nomod' (i.e., background) measurement is given, use it to get rid of strong background spots.
Z_NOMOD_CUTOFF = 6.0;
if ~isempty( d_nomod )  

  b = d_nomod.area_peak;
  if ( size( b, 1 ) ~= size( a, 1 ) )
    fprintf( 'Hey, d_nomod needs to have the same dimensions as d !!!\n' ); return;
  end
  b = quick_norm( b, normbins );

  for i = 1:size( b, 1 )
    Zscore_nomod(i,:) =  ( b(i,:) - mean( b(i,:)) )/ std( b(i,:),0,2);
  end  

  for i = 1:size(b,1)
    for j = 1:size(b,2)
      if Zscore_nomod(i,j) > Z_NOMOD_CUTOFF
	Zscore( i,j ) = 0.0;
      end
    end
  end

end


%% Need to reshape into NRES by NRES matrix.
% Convert factor
ZSCORE_SCALING = -1.0;

NRES = length(d.sequence);
Zscore_full = zeros( NRES, NRES );
pos1 = d.seqpos - d.offset;
pos2 = d.mutpos - d.offset;

% look for first library, as marked by MutPos. Note that first column is assuemd to be wt.
for i = 2 : length( d.mutpos )
  if (i < length( d.mutpos )) & ( isnan( d.mutpos(i+1) ) | d.mutpos(i) > d.mutpos(i+1) ); 
    break; 
  end
end
gp = [2:i];

Zscore_full( pos1, pos2(gp) ) = Zscore(:,gp) * ZSCORE_SCALING;

mutpos = d.mutpos( gp );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_and_save(Z, seqplot,  outfile, print_stuff );

if ~exist( 'print_stuff' ); print_stuff = 0; end

image( seqplot, seqplot,  (-Z' - 1.5 ) * 64 );
hold on
plot( seqplot, seqplot, 'k' ); 
hold off
colormap( 1 - gray(100) );
h=title( outfile );set(h,'interpreter','none');
xlabel( 'seqpos');
ylabel( 'mutpos' );
save( outfile, '-ascii','Z');

set(gcf,'PaperPositionMode','auto');

if ( print_stuff )
  eps_file = [outfile, '.eps' ];
  fprintf( 'Outputting to postscript file: %s\n', eps_file );
  print( eps_file, '-depsc2', '-tiff' );
end

