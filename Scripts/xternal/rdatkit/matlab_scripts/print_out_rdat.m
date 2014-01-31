function print_out_rdat( r, search_tag );
% print_out_rdat( r, search_tag );
%
%  Input: 
%   r          = name of RDAT file, or RDAT object read into MATLAB.
%   search_tag = [optional] text to look for in file.
%  create a bunch of .eps postscript files for cloud lab RDAT. 
%
% (c) R. Das, Stanford University, 2013.
%
if ischar( r )
 r = read_rdat_file( r );
end

SIGNAL_TO_NOISE_TAG = 0;

data_cols = 1:size( r.reactivity, 2 );
if exist( 'search_tag', 'var' );  data_cols = get_data_cols( r, search_tag ); end;

N = length( data_cols );

set( gcf,'Position',[1          -6         725        1090] );

% number of rows per plot
Nrows = 50;
PAD = 1;

NPLOT = ceil(N / Nrows );
scalefactor = 40 / mean(mean( r.reactivity ) );

set(gcf, 'PaperPositionMode','auto','color','white');

clf; subplot(2,1,1);
set(gca,'position',[0.05 0.35 0.9 0.6])


for i = 1:NPLOT
  
  cla;
  Nmin = ( i - 1 ) * Nrows + 1 - PAD;
  Nmax = ( i ) * Nrows + PAD;

  Nmin = max( Nmin, 1 );
  Nmax = min( Nmax, N );
  
  Nrange = [Nmin:Nmax];
  
  max_pos =  size( r.reactivity, 1 );
  tot_counts = sum( r.reactivity( :, data_cols(Nrange) ), 2 );
  blank_rows = find( tot_counts == 0 );
  if ( ~isempty( blank_rows ) ) max_pos = min( blank_rows )-1; end;
  
  image( scalefactor * r.reactivity( 1:max_pos, data_cols(Nrange) ) ); 
  
  M = length( Nrange );
  set(gca,'xlim', [0.5 M+0.5] );
  
  make_lines( 1:M,'k',0.25 );


  for j = 1:M
    
    % sequences
    for k = 1:length( r.seqpos )
      sequence = r.sequences{ data_cols( Nrange(j) ) };
      idx = r.seqpos(k) - r.offset;
      if idx > 0 & idx <=length( sequence );
	seqchar = sequence( idx );
	h = text( j, k, seqchar,'fontsize',7,'fontweight','bold','color',getcolor( seqchar ),'horizontalalignment','center','clipping','on' );
      end
    end

    set(gca,'ydir','normal','tickdir','out','xtick',[],'fontweight','bold');

    % names.
    anot = r.data_annotations{ data_cols( Nrange(j) ) };
    design_name  = get_tag( anot, 'MAPseq:design_name' );
    project_name = get_tag( anot, 'MAPseq:project' );
    signal_to_noise_tag = '';
    if SIGNAL_TO_NOISE_TAG
      signal_to_noise = get_tag( anot, 'signal_to_noise' );
      signal_to_noise_tag = ['[ signal/noise: ',signal_to_noise,'] '];
    end
    tag = [ signal_to_noise_tag, design_name,' ',project_name ];

    h = text( j, 0, tag );
    set(h,'rotation',90,'horizontalalign','right','fontsize',10,'fontweight','bold','interp','none');

  end
  
  colormap( 1 - gray(100 ))

  text(0, max_pos+2, [num2str(i),'/',num2str(NPLOT)],'verticalalign','bottom','fontweight','bold','fontsize',12 );

  draw_structure_boundaries( Nrange, r, data_cols );
      
  print_suffix = '';
  if exist( 'search_tag','var'); print_suffix = ['_', search_tag]; end;
  print_file_name = sprintf('rdat%s_plot%02d.eps',print_suffix, i );
  print( '-depsc2', print_file_name );
  fprintf( 'Outputted:%s\n ', print_file_name );
  
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function color = getcolor( seqchar )
color = 'k';
switch seqchar
 case 'A'
  color = [1 0.6 0];
 case 'C'
  color = [0 0.5 0];
 case 'G'
  color = [1 0 0];
 case 'U'
  color = [0 0 1];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function   draw_structure_boundaries( Nrange, r , data_cols);

hold on

bounds = [];
plot_structures = {};
structure = '';
for j = 1:length( Nrange )
  structure_j = r.structures{ data_cols( Nrange(j) )};

  if ~strcmp( structure_j, structure )
    structure = structure_j;
    plot_structures = [plot_structures, structure ];
    bounds = [ bounds, j ];
  end
end

%plot_structures = [plot_structures, structure ];
bounds = [bounds, length( data_cols( Nrange) )+1];

start_bound = 1;
for j = 1:length( plot_structures )

  structure = plot_structures{j};
 
  start_bound = bounds(j);
  stop_bound  = bounds(j+1)-1;
  
  [start_bound, stop_bound];
  
  for m = 1:length( r.seqpos )

    q = r.seqpos(m) - r.offset;

    if ( q < 1 ) continue; end;
    if ( q > length( structure ) ) continue; end;
    if ( q == 1 | ( structure(q) == '.' & structure(q-1) ~= '.' ) )
      start_struct = q;   
    end

    if ( q == length(structure) | ( structure(q) == '.' & structure(q+1) ~= '.') )
      stop_struct = q;

      plot( [start_bound-0.5, stop_bound+0.5],  (start_struct-0.5)*[1 1], '-','color',[1 0.7 0],'linewid',2 );
      plot( [start_bound-0.5, stop_bound+0.5],  (stop_struct +0.5)*[1 1], '-','color',[1 0.7 0],'linewid',2 );

      if ( start_bound > 1 | Nrange(1) == 1 )
	plot( [start_bound-0.5, start_bound-0.5],  [start_struct-0.5, stop_struct+0.5], '-','color',[1 0.7 0],'linewid',2 );
      end

      if ( stop_bound < length(Nrange) | Nrange(end) == length( data_cols ) )
	plot( [stop_bound+0.5 , stop_bound+0.5],  [start_struct-0.5, stop_struct+0.5], '-','color',[1 0.7 0],'linewid',2 );
      end
    end
    
  
  end

end


