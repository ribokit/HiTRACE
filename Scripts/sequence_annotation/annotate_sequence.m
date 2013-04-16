function [xsel,seqpos,area_pred] = annotate_sequence( d_align, xsel, sequence_full, ...
					offset, data_types, first_RT_nucleotide, structure );
% ANNOTATE_SEQUENCE - Tool for rapid manual assignment of bands in electropherograms.
%
%  [xsel,seqpos,area_pred] = annotate_sequence( d_align, xsel, sequence_full, offset, data_types, first_RT_nucleotide, structure );
%
%
% Input:
%  d_align        = matrix of aligned electrophoretic traces.
%  xsel           = band positions if you already have them. Give [] if not initialized yet.
%  sequence_full  = sequence of RNA. 
%  offset         = value that is added to sequence index to achieve 'historical'/favorite numbering. [default: 0]
%  data_type      = cell of tags of modification reactions in each trace, e.g., 
%               {'SHAPE','SHAPE','ddTTP'}
% first_RT_nucleotide = integer that gives nucleotide immediately 5' of primer binding site. 
%                         [default: end of sequence]
% structure       = structure in dot/bracket notation [give as '' if unknown] [default: '']
%
% Output:
% xsel      = positions of bands across all lanes.
% seqpos    = sequence numbers that go with each xsel
% area_pred = matrix of zeros and ones that mark band locations for entire assignable sequence. 
%
% (C) R. Das, 2013
%

if nargin == 0;  help( mfilename ); return; end;

% initialize outputs.
if ~exist('xsel');  xsel = []; end
seqpos = [];

if length( xsel ) > 1 & ( xsel(2) > xsel(1) )
  fprintf( 'WARNING! WARNING! WARNING!\n' )
  fprintf( 'You are using the old style of marking, where the bands are marked from 3'' to 5'' stop positions.\n' )
  fprintf( 'This script is switching the order!\n')
  fprintf( 'Outputted seqpos will go from small to large values.\n')
  fprintf( 'Outputted xsel will go from large to small values\n')
  xsel = reverse_sort( xsel );
end

if ~exist('offset');  offset = 0; end
if ~exist('period');  period = 1; end
if ~exist('marks');  marks = []; end
if ~exist('mutpos');  mutpos = []; end
if ~exist('data_types'); data_types = []; end
if ~exist('structure');  structure = ''; end

% 
if exist( 'first_RT_nucleotide' ) | isempty( first_RT_nucleotide )
  sequence = sequence_full( 1 : (first_RT_nucleotide-offset) );
  if length( structure ) > length( sequence ); structure = structure( 1: length(sequence ) ); end;
else
  sequence = sequence_full;
end

% fill out area_pred, based on data_types.
numlanes = size(d_align,2);
if isempty( data_types ) % no input.
  area_pred = zeros( length(sequence), numlanes ); 
elseif iscell( data_types ) % input is a matrix
  for m = length( data_types )+1 : numlanes; data_types{m} = ''; end; % pad to number of lanes. 
  area_pred = get_area_pred( sequence, data_types, offset, structure );
else
  if size( data_types, 1 ) == length( sequence )
    area_pred = data_types;
  else
    fprintf( 'Problem: input area_pred/data_types does not have same size [%d] as sequence [%d]\n', size(data_types,1),  length( sequence ) );  
    return;
  end
  data_types = [];
end

ylim = get(gca,'ylim');
ymin = ylim(1);
ymax = ylim(2);

xlim = get(gca,'xlim');
if ( xlim(1) ~= 0.5 | xlim(2) ~= numlanes+0.5 )
  ymin = 1;
  ymax = size(d_align,1);
  colormap( 1 - gray(100) );
end
axis( [ 0.5 numlanes+0.5 ymin ymax]);

scale_factor = 40/ mean(mean(abs(d_align)));
image( d_align * scale_factor );

set (gcf, 'WindowButtonMotionFcn', @mouseMove);

contrast_factor = 100.0/size(colormap,1);
numlanes = size(d_align,2);

stop_sel = 0;
xsel = reverse_sort( xsel );

annotation_handles = {};

update_plot = 1;
update_ylim = 1;
update_contrast = 1;

set(gcf, 'PaperPositionMode','auto','color','white','pointer','fullcross');

while ~stop_sel

  % try to do 'lazy update'
  if update_plot;     annotation_handles = make_plot( d_align, xsel, sequence, offset, area_pred, annotation_handles );   end
  if update_contrast; colormap( 1 - gray( round( 100/contrast_factor ) ) ); end;
  if (update_ylim | update_plot); do_ylim_update( annotation_handles, ymin, ymax ); end;

  update_plot = 0;
  update_contrast = 0;
  update_ylim = 0;

  title( ['j,l -- zoom. i,k -- up/down. q -- quit. left-click -- select. \newline',...
	  'middle button -- erase.', ' r -- reset.', 'x -- auto-assign. ',...
	  sprintf(' # of Annotation: %d / %d', length(xsel),length(sequence))]);  
  
  %  [yselpick, xselpick, button ]  = ginput(1);

  keydown = waitforbuttonpress;
  button = get( gcf, 'SelectionType' );
  mouse_pos = get( gca, 'CurrentPoint' );
  xselpick = mouse_pos(1,2);
 
  if ( ~keydown ) % mousebutton pressed!
    switch( button  )
     case 'normal' % left-click
      xsel = [xsel xselpick];
      xsel = reverse_sort( xsel );
      update_plot = 1;
     case 'extend' % middle click, or shift-click on Mac
      xsel = remove_pick( xsel, xselpick );
      update_plot = 1;
    end
  else    
    keychar = get(gcf,'CurrentCharacter');

    switch keychar
     case {'p', 'P'}
      [name path] = uiputfile('xsel.mat', 'Write xsel to matlab file');
      if(name)
	save(strcat(path, name), 'xsel' );
      end
     case {'o','O'}
      [name path] = uigetfile('*.mat', 'Read xsel from matlab file');
      if(name)
	a = load(strcat(path, name));
	xsel = [xsel a.xsel];
	xsel = reverse_sort( xsel );
      end
      update_plot = 1;
     case {'q','Q','z','Z'}
      stop_sel = 1;
      if( length(xsel) > 0 && (length(xsel) < length(sequence) || length(xsel) > length(sequence)) )
	%if isstruct(USE_GUI)
	%btn = questdlg(sprintf('Warning: You have assigned only %d band position(s), which is less/more than the total number of bands (%d). How do you want to proceed?', length(xsel), length(sequence)), ...
	%	       'Warning!', 'Return to image', 'Force to go','Force to go');
        %if(~strcmp(btn, 'Force to go'))
	%  stop_sel = 0;
        %end
	%else
	%fprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n')
	%fprintf('Warning: You have assigned only %d band position(s), which is less than the total number of bands (%d).\n', length(xsel), length(sequence));
	%fprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n')
	%end
      end
     case {'e','E'} % doesn't work -- use mouse middle-click to erase.
      xsel = remove_pick( xsel, xselpick );
      update_plot = 1;
     case {'j','J'}
      current_relative_pos =  (xselpick - ymin)/(ymax-ymin);
      yscale = (ymax - ymin)*0.75;
      ymin = xselpick - yscale * (current_relative_pos);
      ymax = xselpick + yscale * ( 1- current_relative_pos);
      update_plot = 1; % for text labels
      update_ylim = 1;
     case {'l','L'}
      current_relative_pos =  (xselpick - ymin)/(ymax-ymin);
      yscale = (ymax - ymin)/0.75;
      ymin = xselpick - yscale * (current_relative_pos);
      ymax = xselpick + yscale * ( 1- current_relative_pos);
      update_plot = 1; % for text labels
      update_ylim = 1;
     case {'i','I'}
      yscale = (ymax - ymin);
      ymin = ymin - yscale*0.05;
      ymax = ymax - yscale*0.05;
      update_plot = 1; % for text labels
      update_ylim = 1;
     case {'k','K'}
      yscale = (ymax - ymin);
      ymin = ymin + yscale*0.05;
      ymax = ymax + yscale*0.05;
      update_plot = 1; % for text labels
      update_ylim = 1;
     case {'b', 'B'}
      yscale = (ymax - ymin);
      ymin = ymin + yscale;
      ymax = ymax + yscale;
      update_plot = 1; % for text labels
      update_ylim = 1;
     case {'t', 'T'}
      yscale = (ymax - ymin);
      ymin = ymin - yscale;
      ymax = ymax - yscale;
      update_plot = 1; % for text labels
      update_ylim = 1;
     case {'1'}
      contrast_factor = contrast_factor * sqrt(2);
      update_contrast = 1;
     case {'2'}
      contrast_factor = contrast_factor / sqrt(2);
      update_contrast = 1;
     case {'r', 'R'}
      xsel = []; % reset
      update_plot = 1;
     case {'x','X'}
      seqpos = length(sequence) - [1:length(xsel)] + 1 + offset;			
      
      % um, a guess. Probably could think of a better one...
      peak_spacing = size( d_align, 1 )/ length(  sequence );
      
      if ~exist( 'area_pred' ) | isempty( area_pred )
	fprintf( 'You need to input data_types or area_pred if you want to use auto-assign!\n' )
      else
	input_bounds = [];
	if length( xsel ) == 2; 
	  input_bounds = sort(xsel); 
	  peak_spacing = ( input_bounds(2) - input_bounds(1) ) / length( sequence );
	end;
	area_pred_reverse = area_pred(end:-1:1,:); % should fix auto_assign to reverse.
	fprintf( 'Running auto-assign. This might take a minute.\n');
	xsel = auto_assign_sequence( d_align, sequence, seqpos, offset, area_pred_reverse, peak_spacing, input_bounds, 0, data_types );
	xsel = reverse_sort( xsel ); % should fix auto_assign to reverse.
      end
      update_plot = 1;
    end        
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length( xsel ) > length( sequence )+1
  fprintf( 'WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!\n');
  fprintf( ' Number of selected positions %d exceeds length of sequence %d!\n', length(xsel),length(sequence) );
  fprintf( 'WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!\n');
end
if length( xsel ) < length( sequence )
  fprintf( 'WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!\n');
  fprintf( ' Number of selected positions %d is less than length of sequence %d!\n', length(xsel),length(sequence) );
  fprintf( 'WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!\n');
end

title '';

seqpos = get_seqpos( sequence, offset, xsel );

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xsel = remove_pick( xsel, xselpick )

if length( xsel ) > 0
  [dummy, closestpick] = min( abs(xsel - xselpick) );
  
  xsel = xsel( [1:(closestpick-1) ...
		(closestpick+1):length(xsel)] ...
	       );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  annotation_handles = make_plot( d_align, xsel, sequence, offset, area_pred, annotation_handles );

% an 'annotation handle' consists of the text, horizontal line, and circle markers
%  that go with each sequence position.
% 
% each handle will have { xsel position, handle for line, handle for text, handle for circle1, ... };
%
% the annotation handles are indexed by their sequence position [without using offset], but
% going backward from the very first one [at length(Sequence) ]. This is to allow
% additional bands that are 5' to the nominal sequence beginning.
%

numlanes = size( d_align, 2 ) ;

seqpos = get_seqpos( sequence, offset, xsel );

% why not just set contrast_factor based on colormap?
% might be much faster.

xlim = get(gca,'xlim');
xscale = abs(xlim(1) - xlim(2));
hold on

checked_handle = zeros( length( annotation_handles ),1 );
for i = length( xsel ):-1:1
  
  seq_idx = seqpos(i) - offset;
  handle_idx = length( sequence ) - seq_idx + 1;
  checked_handle( handle_idx ) = 1;  
  if handle_idx <= length( annotation_handles ) & length( annotation_handles{handle_idx} ) > 0 
    
    % This will move up or down handles that have already been created.
    handles = annotation_handles{ handle_idx };
    
    if ( handles{1} == xsel(i) ) continue; end;
    %fprintf( 'Updating handle: %d\n', handle_idx );
    handles{1} = xsel(i);
    annotation_handles{handle_idx} = handles;
    
    for m = 2:length( handles )
      h = handles{m};
      if isprop( h, 'YData' )
	newposition = get(h,'Ydata');
	newposition = newposition * 0 + xsel(i);
	set( h, 'YData', newposition );
      elseif isprop( h, 'Position' )
	newposition =  get( h, 'Position' );
	newposition( 2 ) = xsel(i);
	set( h, 'Position', newposition );
      end
    end
    
  else % need to create the handle.
              
    %fprintf( 'New handle: %d\n', handle_idx );
    
    handles = {};
    handles{1} = xsel(i);
    
    mycolor = [ 1 0 1]; % magenta for unknown.
        
    hold on
    h = plot( [0.5 numlanes+0.5 ], [xsel(i) xsel(i)], 'color',mycolor ); 
    handles = [ handles, h ];
    
    if ( seq_idx >=1  &  seq_idx <= length(sequence) )

      % eterna colors. :)
      seqchar = sequence( seq_idx   );    
      switch seqchar
       case {'A','a'}
	mycolor = [0 0 1];
       case {'C','c'}
	mycolor = [0 0.5 0];      
       case {'U','T','u','t'}
	mycolor = [1 0.5 0];
       case {'G','g'}
	mycolor = [1 0 0];
      end
      set(h(1),'color',mycolor);
      
      txt_to_show = seqchar;
      seqnum = seqpos(i);
      txt_to_show = [seqchar,num2str( seqnum)];
      
      h = text( 0.5, xsel(i), txt_to_show );        
      set(h,'HorizontalAlignment','right');
      set(h,'fontweight','bold','fontsize',6);%,'clipping','on');
      handles = [handles, h ];

      mark_points = find( area_pred(seq_idx,:) > 0.5 );      
      h = plot( mark_points, mark_points*0 + xsel(i) , 'ro' );
      handles = [handles, h ];
  
    end

    annotation_handles{ handle_idx } = handles;    
  end

end

% delete any handles that are no longer tagged to an 'xsel'
delete_handle_idx = find( ~checked_handle );
for i = 1:length( delete_handle_idx )
  handle_idx = delete_handle_idx(i);
  handles = annotation_handles{ handle_idx };
  for m = 2:length( handles)
    %set( handles{m},'visible','off' );
    delete( handles{m} );
  end
  annotation_handles{ handle_idx } = {};
end

axis off
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = reverse_sort( x );
x = sort( x );
x = x( end:-1:1 );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function seqpos = get_seqpos( sequence, offset, xsel );

%seqpos = length( sequence ) + offset + 1 - [1:length(xsel)];
seqpos = length( sequence ) + offset - length(xsel) + [1:length(xsel)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  do_ylim_update( annotation_handles, ymin, ymax ); 

for i = 1:length( annotation_handles )
  handles = annotation_handles{i};

  m = 3; % text in handle
  if length( handles ) < m; continue; end;
  h = handles{m};

  pos = get( h, 'Position' );
  y = pos(2);
  if ( y >= ymin & y <= ymax )
    set( h, 'visible', 'on' );
  else
    set( h, 'visible', 'off' );
  end
end

ylim( [ymin ymax]);


function mouseMove (object, eventdata)
C = get (gca, 'CurrentPoint');
%title(gca, ['(X,Y) = (', num2str(C(1,1)), ', ',num2str(C(1,2)), ')']);
