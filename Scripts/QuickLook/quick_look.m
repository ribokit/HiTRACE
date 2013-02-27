function [d, da, ylim, labels] = quick_look( dirnames, ylim, trace_subset, signals_and_ref, dye_names, moreOptions )
% QUICK_LOOK  basic read-in script for ABI sequencer files for HiTRACE analysis
%
%  [d, da, ylim, labels] = quick_look( dirnames, ylim, trace_subset, signals_and_ref, dye_names, moreOptions )
%
%   quick_look() runs scripts to read in directories of .ab1 files from ABI 
%    do initial data processing, and give feedback.
%  
%  Required input is:
%
%     dirnames = list of directories names  {'dir1','dir2',...'dirN'}
%
%  The rest of the inputs are optional.
%
%     ylim = vector of two values -- the minimum and maximum time point to show (ymin and ymax) [default is "auto-select"]
%     trace_subset = subset of electropherograms to align and return (e.g., [1:8 12:16])
%     signals_and_ref = which color channels have signal. The last entry is the color channel number with
%                                the reference ladder co-loaded in all samples. (Default: [1 4], i.e. 1 has signal, 4 has reference. )
%     dye_names  = Names of dyes for each color channel, used for leakage correction. Currently we assume 4 channels. Input this as a cell of strings.
%                    (Default: {'FAM','ROX'}). 
%                  Note: the number of dye names should correspond to the number of signals_and_ref. If dye names is given as {} or [],
%                    no leakage correction will be applied.
%     moreOptions= Any of the following to turn off data processing steps: {'noPlotStuff', 'noNormalize', 
%             'noSmoothBaselineSubtract',  'noLeakageCorrection' 'noLocalAlign'}  [Default: run al processing steps]
%
% (C) Rhiju Das, 2009-2011, 2013
%   Thanks to S. Denny & C. Cheng for input in multicolor applications
%

d = []; da = []; labels = {}; 

% make backwards compatible...
if exist( 'trace_subset', 'var' ) & length(trace_subset) == 1 & length( trace_subset ) == 1 & trace_subset > ylim
  % this may be the old-style quick-look
  help( 'quick_look' );
  fprintf( '\nWARNING! WARNING! Are you using the old style quick_look with ymin and ymax as separate arguments?\n')
  fprintf( '  Give them instead in a vector [ymin,ymax].\n Or consider not specifying them -- there  is an auto-range finder now! \n\n\n');
  fprintf( 'Press a key to continue\n' )
  pause
  ylim = [ylim, trace_subset ]; trace_subset = [];  
  if exist( 'signals_and_ref' ) trace_subset = signals_and_ref; end;
  signals_and_ref = [];
end

if ~exist( 'ylim', 'var'); ylim = []; end;
if ~exist( 'signals_and_ref', 'var' ) | isempty( signals_and_ref ); signals_and_ref = [1 4]; end;
if length( signals_and_ref ) < 2; fprintf( 'Must have at least 2 channels specified in signals_and_ref!\n'); return; end;
if ~exist( 'dye_names', 'var' )  | ~iscell( dye_names )
  if length( signals_and_ref ) == 2; 
    dye_names = {'FAM','ROX'};
  elseif length( signals_and_ref ) == 3; 
    dye_names = {'FAM','HEX','ROX'}; 
  elseif length( signals_and_ref ) == 4; 
    dye_names = {'FAM','HEX','TAMRA','ROX'}; 
  end    
end;
if length( dye_names ) > 0 & length( signals_and_ref ) ~= length( dye_names) ; fprintf( 'Length of dye_names must be 0 or match length of signals_and_ref\n'); return; end;

% this creates a cell that has several elements, with blank strings where dye_names were not specified.
% for example, if dyenames is { 'FAM', 'ROX' } and signals_and_ref is [1 4], 
%  we know that there are at least 4 color channels, and that the 1st and 4th are specified so
% dye_names_full = { 'FAM', '', '', 'ROX'} carries that information concisely.
dye_names_full = get_dye_names_full( dye_names, signals_and_ref );

sigcol = signals_and_ref( 1 ); % will generalize this in a bit.
refcol = signals_and_ref( end );
fprintf( 'Assuming reference channel: %d\n', refcol );


% Parse some of these crazy options.
if ~exist( 'moreOptions' ) moreOptions = {}; end;
if ~iscell( moreOptions); moreOptions = { moreOptions }; end;
PLOT_STUFF = 1;
NORMALIZE = 1;
LOCAL_ALIGN = 1;
SMOOTH_BASELINE_SUBTRACT = 1;
for m = 1:length( moreOptions )
  switch moreOptions{m}
   case 'noPlotStuff' 
    PLOT_STUFF = 0; 
   case 'noNormalize' 
    NORMALIZE = 0;
   case 'noSmoothBaselineSubtract'
    SMOOTH_BASELINE_SUBTRACT = 0;
   case 'noLeakageCorrection' 
    dye_names_full = {};
   case 'noLocalAlign' 
    LOCAL_ALIGN = 0;
   otherwise
    fprintf( 'Unrecognized options! %s\n', moreOptions{m} ); 
    return;
  end
end

  
d = [];
da = [];
labels = {};
d_noalign = {};
da_noalign = {};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STAGE 1
% The ABI file format is a special binary file format -- 
% I have written scripts to take care of the readin, and
% to put all the data into one object ("data_all")

if ~iscell( dirnames )
  if isdir( dirnames ) 
    dirnames = { dirnames };
  else 
    % Assume it is a text file with names.
    [dirnames] = textread(dirnames,'%s');
  end
end
tag = dirnames{1};
if tag(end) == '/'; tag = tag(1:end-1); end;

line_pos = [ 0 ];
filepath = '';
[ data_all, filenames_all, data_init, data_length ] = ...
    read_abi_dirs( filepath, dirnames, dye_names_full, PLOT_STUFF );
if ~exist( 'trace_subset' ) | length( trace_subset) == 0;  trace_subset = [ 1 : length( data_all ) ]; end

if length( data_all ) == 0 ; return; end;

line_pos = data_init-1;
numfiles = length( data_all );
labels = filenames_all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STAGE 2
%Quick check... let's plot the raw data.
% The signal we care about is in channel 1.
% A reference ladder (in a different color) shows up in channel 4.

% Copy the appropriate columns to a single matrix, to make plotting easy.
% standard...

%truncate data to the same length (this shouldn't be necessary, but happened once with a mixed-up data set from the PAN facility.
for i = length(trace_subset)
    size_v = size(data_all{i},1);
end
max_size = max(size_v);

d0_signal = [];
d0_reference_ladder = [];


for i = 1:length( trace_subset )  
  if(size(data_all{trace_subset(i)},1) < max_size)
    data_all{trace_subset(i)}(max_size, size(data_all{trace_subset(i)},2)) = 0;
  end
    
  % this is a straightforward subtraction of an offset.
  d0_signal(1:size(data_all{trace_subset(i)}(:,1),1),i)           = baseline_subtract(data_all{trace_subset(i)}(:,sigcol));
  d0_reference_ladder(1:size(data_all{trace_subset(i)}(:,1),1),i) = baseline_subtract(data_all{trace_subset(i)}(:,refcol));
end

AUTOFIND_YLIM = 0;
if isempty( ylim ) 
  AUTOFIND_YLIM = 1;
  [ymin, ymax] = findTimeRange( d0_signal );
  ylim = [ymin,ymax];
end
ymin = ylim(1); ymax = ylim(end);

if PLOT_STUFF
  h = figure(2);
  set(h,'Position',[50,50,600,800]);
  set(gcf, 'PaperPositionMode','auto','color','white');
  clf
  subplot(1,2,1);
  image( 0.05 * d0_signal);
  
  h=title( 'Signal (channel 1)');
  set( h,'interpreter','none' );
  
  axis( [ 0.5 length( trace_subset )+0.5 ymin ymax] );
  set( gca, 'xtick', 1:length( trace_subset ), ...
	    'xticklabel', char( labels{trace_subset}  )  );
  xticklabel_rotate;
  make_dividers( trace_subset, line_pos, ymin, ymax );

  %figure(3)
  subplot(1,2,2);
  image( d0_reference_ladder*2);
  title( 'Reference ladder (channel 4)')
  axis( [ 0.5 length( trace_subset )+0.5 ymin ymax] );
  make_dividers( trace_subset, line_pos, ymin, ymax );
    set( gca, 'xtick', 1:length(trace_subset), ...
	    'xticklabel', char( labels{trace_subset}  )  );
  xticklabel_rotate;
  
  colormap( 1- gray(100));
  print( '-depsc2',[tag,'_Figure2.eps']);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STAGE 3

reflane = 1; % this is a vestige of some old testing stuff.

for i = 1:length(data_init)
  start_index = data_init(i);
  if i < length( data_init )
    final_index = data_init(i+1)-1;
  else
    final_index = length( data_all );
  end
  data_align_group{i} = align_capillaries( ...
      { data_all{[start_index:final_index]} }, refcol, reflane);
end

data_align = align_capillaries_group( data_align_group, refcol, 1, 1);


d = []; da= [];

for i = 1:length( trace_subset); 
  d(:,i)  = baseline_subtract(data_align{trace_subset(i)}(:,1));
  da(:,i) = abs(baseline_subtract(data_align{trace_subset(i)}(:,refcol)));

  if(  ymax > size( d, 1 ) ) 
    fprintf( 'WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!\n')
    fprintf( 'WARNING! You specified a ymax (%5d) greater than trace length (%5d) WARNING!\n', ymax, size(d,1));
    ymax = size( d, 1);
    fprintf( 'WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!\n')
  end

  if NORMALIZE
    %d(:,i) = d(:,i)/mean( data_align{trace_subset(i)}(1500:2500,3) );
    d(:,i) =  d(:,i)/mean( abs(d(ymin:ymax,i)));
    da(:,i) = da(:,i)/mean( abs(da(ymin:ymax,i)));
  else
    d(:,i) = d(:,i);
    da(:,i) = da(:,i);
  end
    
end;

%save d.mat d;
%save da.mat da;

if PLOT_STUFF
  h = figure(3);
  set(h,'Position',[100,100,600,800]);
  set(gcf, 'PaperPositionMode','auto','color','white');
  clf
  
  %subplot(1,2,1);
  image( 50*d );
  axis( [ 0.5 length( trace_subset )+0.5 ymin ymax] );
  set( gca, 'xtick', 1:length( trace_subset ), ...
	    'xticklabel', char( labels{ trace_subset}  )  );
  xticklabel_rotate;

  make_dividers( trace_subset, [], ymin, ymax );
  h=title( [dirnames{1}]);
  set( h,'interpreter','none' )
  colormap(  1 - gray(100) )

  %print( '-depsc2',[tag,'_Figure3.eps']);

  
  %h = figure(4);
  %set(h,'Position',[150,150,600,800]);
  %set(gcf, 'PaperPositionMode','auto','color','white');
  %%subplot(1,2,2);
  %image( 50*da );
  %axis( [ 0.5 length( trace_subset )+0.5 ymin ymax] );
  %set( gca, 'xtick', 1:length( trace_subset ), ...
  %	    'xticklabel', char( labels{ trace_subset }  )  );
  %xticklabel_rotate;
  %axis off
  %title( 'aligned reference ladders');
  %
  %make_dividers( trace_subset, [], ymin, ymax );
  %title( 'Reference ladder (channel 4)')
  %colormap(  1 -gray(100) )
  %figure(2)
  %figure(3)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STAGE 4 --  subtract off a smooth (but
%  not necessarily constant) baseline
figure(4);
if SMOOTH_BASELINE_SUBTRACT;   d = baseline_subtract_v2( d, ymin, ymax);  end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STAGE 5
% align_by_DP -- local refinement through 
% piece-wise linear transformation.
d_before_DP  = d(  [ymin:ymax], : );
da_before_DP = da( [ymin:ymax], : );
if (LOCAL_ALIGN) [d, da] = align_by_DP_using_ref( d_before_DP, da_before_DP ); end;

if PLOT_STUFF
  h = figure(4); clf;
  set(h,'Position',[150,150,600,800]);
  set(gcf, 'PaperPositionMode','auto','color','white');
  clf
  
  %subplot(1,2,1);
  image( 50*d );
  axis( [ 0.5 length( trace_subset )+0.5 1 size(d,1)] );
  %axis( [ 0.5 length( trace_subset )+0.5 ymin ymax] );
  set( gca, 'xtick', 1:length( trace_subset ), ...
	    'xticklabel', char( labels{ trace_subset}  )  );
  xticklabel_rotate;

  make_dividers( trace_subset, [], 1, size( d, 1 ));
  h=title( [dirnames{1}]);
  set( h,'interpreter','none' )
  colormap(  1 - gray(100) )

  print( '-depsc2',[tag,'_Figure4.eps']);

  
  h = figure(5); clf;
  set(h,'Position',[200,200,600,800]);
  set(gcf, 'PaperPositionMode','auto','color','white');
  %subplot(1,2,2);
  image( 50*da );
  axis( [ 0.5 length( trace_subset )+0.5 1 size(d,1)] );
  %xlim( [ 0.5 length( trace_subset )+0.5 ] );
  set( gca, 'xtick', 1:length( trace_subset ), ...
	    'xticklabel', char( labels{ trace_subset }  )  );
  xticklabel_rotate;
  axis off
  title( 'aligned reference ladders');

  make_dividers( trace_subset, [], 1, size( d, 1 ));
  title( 'Reference ladder (channel 4)')
  
  colormap(  1 -gray(100) )
  
  figure(2)
  figure(4)
  

end

if PLOT_STUFF
  fprintf( ['\nCreated: ',tag,'_Figure2.eps\n'] );
  %fprintf( ['Created: ',tag,'_Figure3.eps\n'] );
  fprintf( ['Created: ',tag,'_Figure4.eps\n'] );
end

fprintf( '\n' );
fprintf( 'ymin = %d\n', ymin)
fprintf( 'ymax = %d\n\n', ymax)
if length( dye_names_full ) > 0; fprintf( 'Applied leakage correction for color channels.\n' ); end;
if AUTOFIND_YLIM;                fprintf( 'Used auto-find of ymin, ymax.\n' ); end;
if NORMALIZE;                    fprintf( 'Normalized data based on mean peak intensity.\n' ); end;
if SMOOTH_BASELINE_SUBTRACT;     fprintf( 'Applied subtration of smooth base line.\n' ); end;

labels = labels( trace_subset );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  make_dividers( trace_subset, line_pos, ymin, ymax );

hold on
for i = 1:length( trace_subset );
  plot( 0.5+i*[1 1], [ymin ymax],'k-', 'linew',0.25); 
end
for i = 2:length( line_pos );
  plot( 0.5+line_pos(i)*[1 1], [ymin ymax],'k-','linew',2); 
end
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dye_names_full = get_dye_names_full( dye_names, signals_and_ref );
dye_names_full = {};
if length( dye_names ) == 0; return; end;

%dye_names_full = {'FAM','HEX','TAMRA','ROX'};

for m = 1:length( signals_and_ref )
  dye_names_full{ signals_and_ref(m) } = dye_names{m};
end

