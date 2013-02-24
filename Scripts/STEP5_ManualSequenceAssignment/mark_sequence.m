function [xsel,seqpos,area_pred] = mark_sequence( image_x, xsel, sequence_full, ...
					offset, data_types, primer_binding_site, structure );
% MARK_SEQUENCE - Tool for rapid manual assignment of bands in electropherograms.
%
%  [xsel,seqpos,area_pred] = ( image_x, xsel, sequence_full, offset, data_types, primer_binding_site, structure );
%
% Output:
% xsel      = positions of bands across all lanes.
% seqpos    = sequence numbers that go with each xsel
% area_pred = matrix of zeros and ones that mark band locations for entire assignable sequence. 
%
% Input:
%  image_x  = matrix of aligned electrophoretic traces.
%  sequence = sequence of reverse-transcribed RNA. 
%  offset   = value that is added to sequence index to achieve 'historical'/favorite numbering.
%
% (C) R. Das, 2008-2011, 2013
%
% It should be possible to dramatically accelerate this by not calling make_plot every time.
%

% initialize outputs.
if ~exist('xsel');  xsel = []; end
seqpos = [];

if ~exist('offset');  offset = -999; end
if ~exist('period');  period = 1; end
if ~exist('marks');  marks = []; end
if ~exist('mutpos');  mutpos = []; end
if ~exist('data_types'); data_types = []; end
if ~exist('structure');  structure = ''; end

% 
if exist( 'primer_binding_site' ) | isempty( primer_binding_site )
  sequence = sequence_full( 1 : (primer_binding_site-offset-1) );
else
  sequence = sequence_full;
end

% fill out area_pred, based on data_types.
numlanes = size(image_x,2);
if isempty( data_types ) % no input.
  area_pred = zeros( length(sequence), numlanes ); 
elseif iscell( data_types ) % input is a matrix
  for m = length( data_types )+1 : numlanes; data_types{m} = ''; end; % pad to number of lanes. 
  area_pred = get_area_pred( sequence, data_types, structure );
else
  if size( data_types, 1 ) == length( sequence )
    area_pred = data_types;
  else
    fprintf( 'Problem: input area_pred/data_types does not have same size [%d] as sequence [%d]\n', size(data_types,1),  length( sequence ) );  
    return;
  end
end

colormap( 1 - gray(100) );
ylim = get(gca,'ylim');
ymin = ylim(1);
ymax = ylim(2);
if ( ymin==0 & ymax==1)
  image( image_x );
  ylim = get(gca,'ylim');
  ymin = ylim(1);
  ymax = ylim(2);
end

contrast_factor = 40/ mean(mean(abs(image_x)));

%  Probably should just remove this... not clear if JUST_PLOT_SEQUENCE flag is needed anymore...
%if (JUST_PLOT_SEQUENCE )
%  if isstruct(USE_GUI); axes(USE_GUI.displayComponents); end;
%  make_plot( image_x, xsel, ymin, ymax, sequence, JUST_PLOT_SEQUENCE, ...
%	     contrast_factor, offset, period,marks,mutpos, area_pred);
%  return;
%end

numlanes = size(image_x,2);

stop_sel = 0;
xsel = reverse_sort( xsel );

% not clear if USE_GUI is needed anymore
%if isstruct(USE_GUI)
%  axes(USE_GUI.displayComponents);
%  make_plot( image_x, xsel, ymin, ymax, sequence, JUST_PLOT_SEQUENCE, ...
%	     contrast_factor, offset, period,marks,mutpos, area_pred);
%  uiwait( msgbox( 'Are you ready to interactively annotate the sequence?','Ready?','modal' ) )
%end
while ~stop_sel

  % not clear if USE_GUI is needed anymore
  %if isstruct(USE_GUI); axes(USE_GUI.displayComponents);; end;

  make_plot( image_x, contrast_factor, ymin, ymax, xsel, ...
	     sequence, offset, area_pred );
  
  title( ['j,l -- zoom. i,k -- up/down. q -- quit. left-click -- select. \newline',...
	  'x -- auto-assign. ',...
	  'middle button -- undo.', ' r -- reset.', ' t,b -- page up/down.', sprintf(' # of Annotation: %d / %d', length(xsel),length(sequence))]);
  
  [yselpick, xselpick, button ]  = ginput(1);

 
  if(isempty(button))
      continue;
  end
  switch( button )
   case 1
    xsel = [xsel xselpick];
    xsel = reverse_sort( xsel );
   case 2
    xsel = remove_pick( xsel, xselpick );
  end
  
  switch char(button)

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
   case {'e','E'}
    xsel = remove_pick( xsel, xselpick );
   case {'j','J'}
    current_relative_pos =  (xselpick - ymin)/(ymax-ymin);
    yscale = (ymax - ymin)*0.75;
    ymin = xselpick - yscale * (current_relative_pos);
    ymax = xselpick + yscale * ( 1- current_relative_pos);
   case {'l','L'}
    current_relative_pos =  (xselpick - ymin)/(ymax-ymin);
    yscale = (ymax - ymin)/0.75;
    ymin = xselpick - yscale * (current_relative_pos);
    ymax = xselpick + yscale * ( 1- current_relative_pos);
   case {'i','I'}
    yscale = (ymax - ymin);
    ymin = ymin - yscale*0.05;
    ymax = ymax - yscale*0.05;
   case {'k','K'}
    yscale = (ymax - ymin);
    ymin = ymin + yscale*0.05;
    ymax = ymax + yscale*0.05;
   case {'1'}
    contrast_factor = contrast_factor * sqrt(2);
   case {'2'}
    contrast_factor = contrast_factor / sqrt(2);
   % are a,d,w,s even in use? if not, remove.
   case {'a','A'}
    xsel = 1.005*(xsel-min(xsel))+min(xsel);
   case {'d','D'}
    xsel = (1/1.005)*(xsel-min(xsel))+min(xsel);
   case {'w','W'}
    xsel = xsel - 5;
   case {'s','S'}
    xsel = xsel + 5;
   % these are 'big' jumps. again not in use -- so remove?
   case {'b', 'B'}
    yscale = (ymax - ymin);
    ymin = ymin + yscale;
    ymax = ymax + yscale;
   case {'t', 'T'}
    yscale = (ymax - ymin);
    ymin = ymin - yscale;
    ymax = ymax - yscale;
   case {'r', 'R'}
    xsel = []; % reset
   case {'x','X'}
    seqpos = length(sequence) - [1:length(xsel)] + 1 + offset;			

    % um, a guess. Probably could think of a better one...
    peak_spacing = size( image_x, 1 )/ length(  sequence );
      
    if ~exist( 'area_pred' ) | isempty( area_pred )
      fprintf( 'You need to input data_types or area_pred if you want to use auto-assign!\n' )
    else
      input_bounds = [];
      if length( xsel ) == 2; 
	input_bounds = sort(xsel); 
	peak_spacing = ( input_bounds(2) - input_bounds(1) ) / length( sequence );
      end;
      area_pred_reverse = area_pred(end:-1:1,:); % should fix auto_assign to reverse.
      xsel = auto_assign_sequence( image_x, sequence, seqpos, offset, area_pred_reverse, peak_spacing, input_bounds, 0 );
      xsel = reverse_sort( xsel ); % should fix auto_assign to reverse.
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
function  make_plot( image_x, contrast_factor, ymin, ymax, ...
		     xsel, sequence, offset, area_pred );

numlanes = size( image_x, 2 ) ;

seqpos = get_seqpos( sequence, offset, xsel );

% why not just set contrast_factor based on colormap?
% might be much faster.
image(contrast_factor * abs(image_x));

xlim = get(gca,'xlim');
xscale = abs(xlim(1) - xlim(2));
axis( [ 0.5 numlanes+0.5 ymin ymax]);

hold on
xsel_to_plot = find( xsel >= ymin & xsel <= ymax );

for i = xsel_to_plot

  mycolor = [ 1 0 1]; % magenta for unknown.

  h1 = plot( [0.5 numlanes+0.5 ], [xsel(i) xsel(i)], 'color',mycolor ); 
  
  seq_idx = seqpos(i) - offset;
  
  if seq_idx < 0 |  seq_idx > length(sequence); continue; end;
  
  seqchar = sequence( seq_idx   );    
  txt_to_show = seqchar;
  seqnum = seqpos(i);
  txt_to_show = [seqchar,num2str( seqnum)];
  
  h2 = text( 0.5, xsel(i), txt_to_show );        
  set(h2,'HorizontalAlignment','right');
  
  for j = 1:numlanes
    if ( area_pred(seq_idx,j) == 1)	    
      plot( j, xsel(i), 'ro' );
    end
  end
  
  % eterna colors. :)
  switch seqchar
   case {'A','a'}
    mycolor = [0 0 1];
   case {'C','c'}
    mycolor = [0 0.7 0];      
   case {'U','T','u','t'}
    mycolor = [1 0.5 0];
   case {'G','g'}
    mycolor = [1 0 0];
  end
  
  set(h1,'color',mycolor);
  set(h2,'fontweight','bold','fontsize',6,'clipping','off');
  

end
  
axis off
hold off;

axis( [ 0.5 numlanes+0.5 ymin ymax]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = reverse_sort( x );
x = sort( x );
x = x( end:-1:1 );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function seqpos = get_seqpos( sequence, offset, xsel );

%seqpos = length( sequence ) + offset + 1 - [1:length(xsel)];
seqpos = length( sequence ) + offset - length(xsel) + [1:length(xsel)];
