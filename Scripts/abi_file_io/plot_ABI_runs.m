function [ data, filenames ] = plot_ABI_runs( dirname, dye_names_full, PLOT_STUFF )
% PLOT_ABI_RUNS: Read in .ab1 files from a directory, correct for cross-channel contamination.
%
%  [ data, filenames ] = plot_ABI_runs( dirname, dye_names_full, PLOT_STUFF )
%
% INPUTS:
%  dirname    = directory with ABI files [.ab1 or .fsa format]
%  dye_names  = [optional] names of dyes in each color channel. default = {}, which means no leakage correction.
%                   Can also specify filename with leakage matrix.
%  PLOT_STUFF = [optional, ignore for now] default: 1.
%
% OUTPUTS:
%  data      = cell containing data matrices for all .ab1/.fsa files. One matrix for each color channel.
%  filenames = which .ab1/.fsa files were read in.
%
% (C) R. Das, 2008-2011, 2013
%

data = {};
filenames = {};
if nargin == 0;  help( mfilename ); return; end;

if ~exist( 'dye_names' ) dye_names = {}; end;
if ~exist( 'PLOT_STUFF' );  PLOT_STUFF = 1; end
if ~exist('ymax')
  ymax=2000;
end

datafiles = dir( [dirname,'/*ab1'] );
if length( datafiles ) == 0
% fragment analysis files -- appear to be the same format as ab1
  datafiles = dir( [dirname,'/*fsa'] ); 
end

count = 0;

data_in = {};
for k = 1:length( datafiles );
  count= count + 1;
  datafile = datafiles( k ).name;
  filenames_in{ count } = datafile( 1:(end-4) ); %remove .ab1 tag.
  datafile = [dirname,'/',datafile];
  filenames_full{ count } = datafile;
  d = read_abi( datafile );  
  data_in{ count } =  d;
end

%clf; plot( data_in{1}(:,1),'r' ); pause;

if length( data_in ) == 0; 
  fprintf( 'WARNING!!! Could not read .ab1 files from: %s\n', dirname );
  return; 
end; % did the files exist?


% in case there are more color channels than specified in dye_names_full...
if ~ischar( dye_names_full )
  for i = length( dye_names_full)+1 : size( d,2 )
    dye_names_full{ i } = '';
  end
end

% leakage correction. Will not do anything if dye_names wasn't specified.
if ischar( dye_names_full )
  lm = load( dye_names_full );
else
  lm = get_leakage_matrix( dye_names_full );
end
data_correct = correct_leakage( data_in, lm );
data_in = data_correct;

%how_many_cols = figure_out_cols( filenames );
for k = 1:count
  whichwell(k) = figure_out_well(filenames_full{k} );
end

[whichwell_sort, reindex ] = sort( whichwell );

for k = 1:count
  data{k} = data_in{reindex(k)};
  filenames{k} = filenames_in{reindex(k)};
end

colorcode = [0 0 1; 0 0.5 0; 1 0.5 0; 1 0 0];


if PLOT_STUFF
  h = figure(1);
  set(h,'Name','Raw traces');
  set(h,'Position',[10,10,600,800]);
  set(gcf, 'PaperPositionMode','auto','color','white');
  clf;
  for k = 1:count
    
    whichwell_mod16 = mod(whichwell_sort(k) - 1,16) + 1;
    
    whichrow = mod(whichwell_mod16-1,8)+1;
    whichcol = floor( (whichwell_mod16-1) /8) + 1;
    whichplot = (whichrow-1)*2 + whichcol;
    
    subplot(8,2,whichplot );
    cla
    xmin = 0; xmax = 8500;
    ymin = -0.5*ymax;
    
    plot( [xmin xmax],[0 0],'k');
    d = data{k};
    %d(:,3) = d(:,3) * 3;
    for i = 1:4
      plot( d(:,i),'color',colorcode(i,:) );
      hold on
    end
    
    h=text( xmin,ymax,filenames{ k },'fontweight','bold' );
    set(h,'interpreter','none','fontsize',8);
    axis( [xmin xmax ymin ymax]);
    
    %ymax = median(sum(data{i},2));
    set(gca,'xtick',[],'box','off')
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ how_many_cols ] = figure_out_cols( filenames );

how_many_cols = 1;
letters = 'ABCDEFGH';

for i = 1:length(letters)
  count = 0;
  
%  searchstring = ['/',letters(i)];
  searchstring = ['\',letters(i)]; % DOS

  for j = 1:length( filenames )
    findstring = strfind( filenames{j}, searchstring);
    if ( ~isempty( findstring ) )
      for k = 1:length( findstring )
	found_well = 0;
	if (filenames{j}(findstring(k)+2) =='0' | ...
	    filenames{j}(findstring(k)+2) == '1' )
	  found_well = 1;
	end
	if found_well
	  count = count+1;
	end
	%colnum = filenames{j}( findstring+2 );
	%whichcols{count} = colnum ;  
      end
    end
  end
  if ( count > 1 )
    how_many_cols = 2;
    return;
  end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  whichplot = figure_out_subplot(filename, how_many_cols );

letters = 'ABCDEFGH';

whichplot = 1;

for i = 1:length( letters )
  %searchstring = ['/',letters(i)];
  searchstring = ['\',letters(i)]; %DOS
  findstring = strfind( filename, searchstring);
  if ( ~isempty( findstring ) ) 
    if ( how_many_cols == 1 )
      for k = 1:length( findstring )
	if (filename(findstring(k)+2) =='0' | ...
	    filename(findstring(k)+2) == '1' )
	  whichplot = i;
	  return;
	end
      end
    else      
      for k = 1:length( findstring )
	if (filename(findstring(k)+2) =='0' | ...
	    filename(findstring(k)+2) == '1' )
	  colnum = filename( findstring(k)+3 );
	  col =  mod( str2num(colnum) - 1, 2 ) + 1;
	  whichplot = 2 * (i-1) + col;
	  return;
	end
      end
    end
  end
end



return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  whichplot = figure_out_well(filename )

whichplot = figure_out_well_for_char_before_searchstring( filename, '/' );
if ( whichplot > 0 ); return; end;

whichplot = figure_out_well_for_char_before_searchstring( filename, '\' );
if ( whichplot > 0 ); return; end;

whichplot = figure_out_well_for_char_before_searchstring( filename, '_' );
if ( whichplot > 0 ); return; end;

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function whichplot = figure_out_well_for_char_before_searchstring( filename, char_before_searchstring );

letters = 'ABCDEFGH';

whichplot = 0;
% Works for ABI3100
for i = 1:length( letters )
  searchstring = [char_before_searchstring,letters(i)];
  %searchstring = ['\',letters(i)]; %DOS
  findstring = strfind( filename, searchstring);
  if ( ~isempty( findstring ) ) 
    for k = 1:length( findstring )
      if ((filename(findstring(k)+2) =='0' | ...
        filename(findstring(k)+2) == '1' )&& ...
        ~isempty(str2num(filename(findstring(k)+3)))&& ...
        isempty(str2num(filename(findstring(k)+4))))
    	colnum = filename( findstring(k)+2:findstring(k)+3 );
        %col =  mod( str2num(colnum) - 1, 2 ) + 1;
        col =  str2num(colnum);
        whichplot = i + 8 *(col-1);
        return;
      end
    end
  end
end
				  
return