function leakage_matrix = figure_out_leakage_correction( infile, outfile, ylimit )
%
% leakage_matrix = figure_out_leakage_correction( dirname, outfile, ylimit )
%
%   OR [to create one-by-one text files]
%
% leakage_vector = figure_out_leakage_correction( infile, outfile, ylimit )
%
% (C) R. Das, 2013.

if nargin == 0;  help( mfilename ); return; end;

leakage_matrix = [];

if ~exist( 'ylimit', 'var' ) ylimit = []; end;% will try to figure this out.

if ischar( infile ) & length( infile ) > 4 &  ...
      ( strcmp( infile(end-3:end), '.ab1') | strcmp( infile(end-3:end), '.fsa') )
  % just one file
  
  leakage_vector = figure_out_leakage_correction_one_file(  infile, outfile, ylimit );
  leakage_matrix = leakage_vector;
  
else   % directory of files.
  
  % go through each file.
  dirname = infile;
  datafiles = dir( [dirname,'/*ab1'] );
  if length( datafiles ) == 0;    datafiles = dir( [dirname,'/*fsa'] );   end

  for k = 1:length( datafiles )
    infile = [dirname, '/', datafiles( k ).name];
    leakage_matrix(:,k) = figure_out_leakage_correction_one_file(  infile, '', ylimit );
    pause;
  end

  if size( leakage_matrix,1) ~= size( leakage_matrix, 2)
    fprintf( 'WARNING!! Number of .ab1/.fsa files (%d) does not match number of channels (%d)!\n', size( leakage_matrix,2), size( leakage_matrix,1) );
  end
  
  fid = fopen( outfile, 'w' );
  for k = 1:length( datafiles )
    fprintf( fid, '\t%8.3f', leakage_matrix(:,k) );
    fprintf( fid, '\n');
  end
  fclose( fid );

  fprintf( 'Outputted: %s\n', outfile );
  
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  leakage_vector = figure_out_leakage_correction_one_file(  infile, outfile, ylimit );

leakage_vector = [];

d = read_abi( infile );
N = size( d, 2 );

% subtract a linear offset 
d_sub = baseline_subtract( d );

if ~exist( 'ylimit', 'var' ) | isempty( ylimit) % try to figure this out.
  [dummy, ypeak] = max( sum( d_sub' ) );
  ymin = ypeak - 50;
  ymax = ypeak + 50;
  ymin = max(ymin,1);
  ymax = min(ymax, size(d,1) );
else
  ymin = ylimit(1);
  ymax = ylimit(2);
end
d_fit = d_sub( ymin:ymax,:);

[dummy, max_channel] = max( sum( d_fit ) );
L = size( d_fit, 1 );

basis_set = [ d_fit(:,max_channel), ones(L,1) ];
for i = 1:N;
  X = lsqr( basis_set, d_fit(:,i) );
  subplot(N,1, i)
  plot( d_sub(ymin:ymax,i),'ko'); hold on
  plot( basis_set*X,'r' ); hold off  
  leakage_vector(i) = X(1);
  h = title( sprintf('%s: Channel %d',infile, i) );
  set(h,'interp','none');
end
legend( 'data',['fit using channel %d',max_channel]);

leakage_vector = leakage_vector/ sum( leakage_vector );

if length( outfile ) > 0
  fid = fopen( outfile, 'w' );
  fprintf( fid, '\t%8.3f', leakage_vector);
  fprintf( fid, '\n');
  fclose( fid );
  
  fprintf( 'Outputted: %s\n', outfile );
end

set(gcf, 'PaperPositionMode','auto','color','white');



