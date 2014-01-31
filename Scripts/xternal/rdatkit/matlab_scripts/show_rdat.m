function rdat = show_rdat(rdat, print_postscript );
%  show_rdat(rdat )
%
%   rdat is 'rdat' file format for high-throughput RNA capillary electrophoresis data
%   label_step  puts a label marker at positions every label_step. (default is 1)
%   show_xsels  shows where peak positions were marked on the 'raw' electrophoretic traces. (default is 0)
%
% Copyright P. Cordero, R. Das, Stanford University, 2010.
%

if nargin==0; help( mfilename ); return; end;

filename = '';
if ischar( rdat )
  filename = rdat;
  rdat = read_rdat_file( filename );
end
if ~exist( 'label_step' );  label_step = 1; end
if ~exist( 'show_xsels' ); show_xsels = 1; end
if ~exist( 'print_postscript' ); print_postscript = 0; end

d = rdat.reactivity;
t = rdat.trace;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up some labels.
xlab = {};
xlab_full = {};
ylab = {};
ylabpos = rdat.seqpos - rdat.offset;
for j = 1:length( ylabpos )
    ylab = [ylab, strcat( num2str( rdat.seqpos(j) ), rdat.sequence( ylabpos(j) ) )];
end 

for j=1:size(rdat.reactivity,2);  xlab{j} = ''; end

MAX_LABEL_LEN = 40;
if length( rdat.data_annotations ) > 0
  for j=1:length(rdat.data_annotations )
    d_annot = rdat.data_annotations{j};
    xlab{j} = remove_tag_before_colon_cell( d_annot );
    xlab_full{j} = strrep(xlab{j},' ',', ');
    % truncate? 
    max_x =  min( length( xlab{j} ), MAX_LABEL_LEN );
    xlab{j} = xlab{j}(1:max_x);
    for i = (length( xlab{j} ) +1) : MAX_LABEL_LEN; xlab{j} = [' ', xlab{j} ]; end;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% show values, i..e. fitted peak intensities or reactivity.
h = figure(1); clf;
set(h, 'Name', 'RDAT Preview');
set(h, 'Position', [100, 100, 800, 600]);
set(h, 'PaperOrientation', 'Landscape', 'PaperPositionMode', 'Manual', ...
      'PaperSize', [11 8.5], 'PaperPosition', [-0.65 0.15 12 8], 'Color', 'White');
colormap( 1 - gray(100));

d_filter = filter_ERROR_lanes( d, rdat.data_annotations );
seq_order = [1:length( rdat.seqpos )];
if ( length(rdat.seqpos) > 1 & rdat.seqpos( 2 ) < rdat.seqpos( 1 ) ) seq_order = length( rdat.seqpos ):-1:1;end;

image( 40 * d_filter( seq_order, :)' /mean(mean(max(d_filter,0)))  )

% sorry for the x-y switch; decided in the end to transpose everything.
if length(ylabpos) > 48; x_ft_sz = 6; else; x_ft_sz = 10; end;
set(gca, 'fontsize', x_ft_sz);
set(gca, 'xTick', [1:length(ylabpos)],'xTickLabel', char(ylab(seq_order)) );
xticklabel_rotate_rdat();
% after xticklabel_rotate, can use xticklabels again to draw structure
%if length(rdat.structure) > 0
%  set(gca,'xTick',[1:length(ylabpos)],'xTickLabel', rdat.structure( ylabpos(seq_order) )','xaxislocation','top' );
%end

if size(xlab,2) > 48; y_ft_sz = 6; else; y_ft_sz = 10; end;
set(gca, 'yTick', [1:size(xlab,2)], 'yTicklabel', char(xlab),'yaxisloc','left');
plot_title = [ rdat.name,'; ', remove_tag_before_colon_cell( rdat.annotations )];
if length( filename ) > 0; plot_title = [ filename,': ',plot_title ]; end;
h = title( plot_title );
set(gca,'fontsize',y_ft_sz,'fontweight','bold','tickdir','out','ticklength',[0 0]);
set( h, 'interpreter','none','fontsize',10,'fontweight','bold');

if length( rdat.structure > 0 ) & length( strfind( rdat.structure, '(') ) > 0 & length( rdat.seqpos ) > 1
  hold on
  for i = 1:length( seq_order )
    if rdat.structure( rdat.seqpos( seq_order(i) ) - rdat.offset ) == '.'
      plot( [i i], [0, size( d_filter, 2 )+1], 'k', 'linewidth', 0.25 );
    end
  end
  hold off
end

if size( d_filter, 2) < 200; make_lines_horizontal( [0:1:size(d_filter,2)],'k',0.25  ); end;
  
if ( print_postscript & length( filename )  > 0 ); 
  eps_file = [filename,'.eps']; fprintf( 'Outputting: %s\n',eps_file );
  print( eps_file, '-depsc2','-tiff' );
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% show aligned 'traces' or electropherograms.
if length( t ) > 0
  figure(2)
  set(gcf,'color','white');
  clf;  colormap( 1 - gray(100));

  t_filter = filter_ERROR_lanes( t, rdat.data_annotations );
  image( 100 * t_filter'/mean(mean(t_filter)))
  hold on;

  if ( length( rdat.xsel_refine ) > 0  | length( rdat.xsel ) > 0 )
    if length( rdat.xsel_refine ) > 0  & show_xsels
      x = rdat.xsel_refine;
      for i=1:size(x,1)      
	plot( [x(i,1), x(i,:), x(i,end)],...
	      [0.5, 1:size(x,2), size(x,2)+0.5 ], 'b'  ); 
	plot( x(i,:),...
	      [1:size(x,2)], '.','color','b'  ); 
      end
      x = rdat.xsel_refine( :, end );
    elseif length( rdat.xsel > 0 )
      x = rdat.xsel;
    end
    [ytickpos, reorder] = sort(x);
    if length(ytickpos) > 48; y_ft_sz = 6; else; y_ft_sz = 10; end;
    if length(xlab) > 48; x_ft_sz = 6; else; x_ft_sz = 10; end;

    set(gca,'fontsize',y_ft_sz);
    set(gca, 'xTick', ytickpos, 'xTickLabel', char(ylab(reorder)),'fontweight','bold');
    set(gca,'fontsize',x_ft_sz);
    set(gca, 'yTick',[1:length(xlab)],'yTicklabel',char(xlab));
    ticklength = get(gca,'ticklength');
    ticklength = [ 0.0 0.0 ];
    set(gca,'ticklength',ticklength,'tickdir','out' );
    
    set( gca,'xdir','reverse');

    xticklabel_rotate_rdat();
  end

  make_lines_horizontal( [0:1:size(t,2)],'k',0.25  )

  %title( ['Traces: ', rdat.name ] );

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% show bargraph with errors
num_lanes = size( rdat.reactivity, 2 );
if ~isempty( rdat.reactivity_error ) && num_lanes <= 16 % totally arbitrary cutoff;
  h = figure(3);
  set(h, 'Name', 'RDAT Plot');
  set(h, 'Position', [100, 100, 800, 600]);
  set(h, 'PaperOrientation', 'Landscape', 'PaperPositionMode', 'Manual', ...
      'PaperSize', [11 8.5], 'PaperPosition', [-0.65 0.15 12 8], 'Color', 'White');
    
  for i = 1: num_lanes
    subplot( num_lanes, 1, i );
        
    % this is a bit involved -- trying to set axes, but ignoring outliers.
    vals = sort( max(rdat.reactivity(:,i),0) );
    minval = vals( round(length(vals)/20) ); % 5th percentile
    maxval = vals( round(length(vals)*19/20) ); % 95th percentile
    valrange = maxval - minval;
    maxval = maxval + 0.6*valrange;
    minval = minval - 0.1*valrange;

    min_x = min(rdat.seqpos) - 1;
    max_x = max(rdat.seqpos) + 1;

    %errorbar( rdat.seqpos, rdat.reactivity(:,i), rdat.reactivity_error(:,i) )
    bar( rdat.seqpos, rdat.reactivity(:,i), 'facecolor',[0.5 0.5 0.5] );%, rdat.reactivity_error(:,i) )
    hold on
    for j = 1:length( rdat.seqpos )
      plot( rdat.seqpos(j) + [0 0 ], rdat.reactivity(j,i) + rdat.reactivity_error(j,i)*[-1 1], 'k','linewidth',2 )
    end
    
    if ( length( rdat.structure ) > 0 )
      structure_plot = zeros( 1, length( rdat.structure ) );
      structure_plot( strfind( rdat.structure, '.' ) ) = 0.7 * maxval;
      plot( [1:length(rdat.structure)]+rdat.offset, structure_plot, '-','color',[0.6 0.3 0.3] );
    end

    plot( [min_x max_x ], [0 0 ],'k' );
    hold off
   
    set(gca,'fontsize',12,'fontweight','bold');
    if (maxval > minval )
      axis( [ min_x max_x minval maxval ] );
    end
    colormap( 1 - gray(100) );
    %if length( xlab ) >= i; title( xlab{i} ); end;

    if length( xlab_full ) >= i;
      h = text( min_x, maxval+(maxval-minval)*0.2, xlab_full{i} ); 
      set(h,'interpreter','none','fontsize',8);
    end
    %set(gca,'xgrid','on','ytick',[]);
    axis off
  end

  axis on; box off
  set(gca,'ytick',[]);
  %xlabel( 'sequence position' );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stripping off the 'mutation' tag makes the labels easier to read.
% If there is more than one annotation, we should fix this function...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tag = remove_tag_before_colon( annotations )
[t,r] = strtok( annotations, ':' );
if strcmp(t,'EteRNA') | strcmp(t,'MAPseq' )
  [t,r] = strtok( r, ':' );
end
%if strcmp(t,'tag')
%  [t,r] = strtok( t, ':' );
%end
tag = r(2:end);

%colon_pos = strfind( annotations, ':' );
%if ~isempty( colon_pos )
%  tag = annotations((colon_pos(1)+1):end);
%end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tag = remove_tag_before_colon_cell( annotations )
tag = '';
for i = 1:length( annotations )
  if ( length( annotations{i} ) > 8 & strcmp( annotations{i}(1:8), 'sequence' ) ) ; continue; end;
  if ( length( annotations{i} ) > 9 & strcmp( annotations{i}(1:9), 'structure' ) ); continue; end;

  tag = [ tag, ' ',  remove_tag_before_colon( annotations{i}) ];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d_filter = filter_ERROR_lanes( d, data_annotations );

d_filter = d;

if length( data_annotations ) < 1; return; end
if length( data_annotations ) ~= size(d,2); fprintf( 'Warning! Mismatch involving data_annotations!\n' ); return; end;

for i = 1:size( d, 2 )
  ok = 1;
  for j = 1:length( data_annotations{i} )
    if strfind( data_annotations{i}{j}, 'ERROR' );
      ok = 0; break;
    end
    if strfind( data_annotations{i}{j}, 'badQuality' );
      ok = 0; break;
    end
  end
  if ~ok
    d_filter(:,i) = d_filter(:,1); 
  end;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = cell2str( c, delim )
if  ~exist('delim','var'); delim = ' '; end;
s = '';
if length(c) > 0
  s = c{1};
  for i = 2:length(c); s = [s,delim,c{i}]; end
end
