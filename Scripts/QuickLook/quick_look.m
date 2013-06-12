function [d, d_ref, ylimit, labels] = quick_look( dirnames, ylimit, trace_subset, signals_and_ref, dye_names, lane_names, moreOptions )
% QUICK_LOOK  basic read-in script for ABI sequencer files for HiTRACE analysis
%
%  [d, d_ref, ylimit, labels] = QUICK_LOOK( dirnames, ylimit, trace_subset, signals_and_ref, dye_names, lane_names, moreOptions )
%
%   quick_look() runs scripts to read in directories of .ab1 files from ABI 
%    do initial data processing, and give feedback.
%  
%  Required input is:
%
%     dirnames = list of directories names  {'dir1','dir2',...'dirN'}
%
%  The rest of the inputs are optional -- specify [] if you want the default.
%
%     ylimit = vector of two values -- the minimum and maximum time point to show (ymin and ymax) [default is "auto-select"]
%     trace_subset = subset of electropherograms to align and return (e.g., [1:8 12:16])
%     signals_and_ref = which color channels have signal. The last entry is the color channel number with
%                                the reference ladder co-loaded in all samples. (Default: [1 4], i.e. 1 has signal, 4 has reference. )
%     dye_names  = Names of dyes for each color channel. Input this as a cell of strings.  Example: {'FAM','ROX'}). 
%                  Note: the number of dye names should correspond to the number of signals_and_ref. If dye names is given as {},
%                    no leakage correction will be applied. (Default: no leakage correction.)
%                  You can also input your own leakage matrix here (give the filename as a string).
%     lane_names = Names of lanes to show on x-axis of each figure. Input this as a cell of strings.
%                  Note: the number of lane names should correspond to the number of total lanes in trace_subset.
%                    If lane names is given as {}, labels from ABI files will be used.
%     moreOptions= Any of the following to turn off data processing steps: {'noNormalize', 
%             'noSmoothBaselineSubtract',  'noLeakageCorrection' 'noLocalAlign'}  [Default: run al processing steps]
%
%  Outputs:
%
%     d      = all signal traces, after baseline subtraction, normalization & alignment. If user asks for more than one channel, 
%               traces for first channel requested are given first, then traces for second channel, etc.
%     d_ref  = reference trace, after all of above. Should be aligned ladders.
%     ylimit = [ymin,ymax] which were found by auto-selection (or input by user)
%     labels = plot labels for each trace, extracted from ABI files.
%
%
% (C) Rhiju Das, 2009-2011, 2013
% minor modified by CYC, T47, May 2013.
%

if nargin == 0;  help( mfilename ); return; end;
d = []; d_ref = []; labels = {}; 

% make backwards compatible...
if exist( 'trace_subset', 'var' ) && length(trace_subset) == 1 && length( trace_subset ) == 1 && trace_subset > ylimit
  % this may be the old-style quick-look
  help( 'quick_look' );
  fprintf( '\nWARNING! WARNING! Are you using the old style quick_look with ymin and ymax as separate arguments?\n')
  fprintf( '  Give them instead in a vector [ymin,ymax].\n Or consider not specifying them -- there  is an auto-range finder now! \n\n\n');
  fprintf( 'Press a key to continue\n' )
  pause
  ylimit = [ylimit, trace_subset ]; trace_subset = [];  
  if exist( 'signals_and_ref' ) trace_subset = signals_and_ref; end;
  signals_and_ref = [];
end

if ~exist( 'ylimit', 'var'); ylimit = []; end;
if ~exist( 'signals_and_ref', 'var' ) || isempty( signals_and_ref ); signals_and_ref = [1 4]; end;
if length( signals_and_ref ) < 2; fprintf( 'Must have at least 2 channels specified in signals_and_ref!\n'); return; end;
if (~exist( 'dye_names', 'var' )  || isempty( dye_names )) 
  dye_names = {}; % signal to not apply a leakage correction -- later use FAM/ROX as default.
end
%if ~exist( 'dye_names', 'var' )  | ~iscell( dye_names )
%  if length( signals_and_ref ) == 2; 
%    dye_names = {'FAM','ROX'};
%  elseif length( signals_and_ref ) == 3; 
%    dye_names = {'FAM','HEX','ROX'}; 
%  elseif length( signals_and_ref ) == 4; 
%    dye_names = {'FAM','HEX','TAMRA','ROX'}; 
%  end    
%end;
if ~ischar( dye_names ) && ~isempty( dye_names ) && length( signals_and_ref ) ~= length( dye_names) ; fprintf( 'Length of dye_names must be 0 or match length of signals_and_ref\n'); return; end;

% this creates a cell that has several elements, with blank strings where dye_names were not specified.
% for example, if dyenames is { 'FAM', 'ROX' } and signals_and_ref is [1 4], 
%  we know that there are at least 4 color channels, and that the 1st and 4th are specified so
% dye_names_full = { 'FAM', '', '', 'ROX'} carries that information concisely.
dye_names_full = get_dye_names_full( dye_names, signals_and_ref );
sigchannels = signals_and_ref( 1 : end-1 ); % will generalize this in a bit.
refchannel = signals_and_ref( end );
fprintf( 'Assuming reference channel: %d\n', refchannel );

% format lane_names to correct length
if ~exist('lane_names','var') || isempty(lane_names); 
    lane_names = {}; 
elseif length(lane_names) < length(trace_subset);
    lane_l_org = length(lane_names);
    lane_names= [lane_names cell(1, length(trace_subset) - lane_l_org)];
    lane_names((lane_l_org + 1):length(trace_subset)) = {' '};
elseif length(lane_names) > length(trace_subset);
    lane_names = lane_names(1:length(trace_subset));
end;

% Parse some of these crazy options.
if ~exist( 'moreOptions','var' ); moreOptions = {}; end;
if ~iscell( moreOptions); moreOptions = { moreOptions }; end;
FIX_LEAKAGE_OF_SATURATING_SIGNALS_TO_REF = 1;
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
d_ref = [];
labels = {};
d_noalign = {};
d_ref_noalign = {};

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
if tag(end) == '/'; tag = tag(1:end-1); end; % get rid of final slash
%  just get the last directory name.
tags = split_string( tag, '/' );
tag = tags{ end };

filepath = '';
[ data_all, filenames_all, data_set_starts, data_length ] = ...
    read_abi_dirs( filepath, dirnames, dye_names_full, PLOT_STUFF );        %Reads in ABI files, performs leakage correction if specified, then plots traces

if isempty( data_all ) ; return; end;

subset_pos = data_set_starts - 1;
numfiles = length( data_all );
filenames_all;

if ~exist( 'trace_subset','var' ) || isempty( trace_subset);  
  trace_subset = 1 : length( data_all );
else
  subset_pos = 0; % if user asks for a subset, it will be hard to define boundaries between subsets.
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STAGE 2
%Quick check... let's plot the raw data.
% The signal we care about is in channel 1.
% A reference ladder (in a different color) shows up in channel 4.

% Copy the appropriate columns to a single matrix, to make plotting easy.
% standard...

%truncate data to the same length (this shouldn't be necessary, but happened once with a mixed-up data set from the PAN facility.
for i = 1:length(trace_subset)
  size_v(i) = size( data_all{trace_subset(i)}, 1);
end
max_size = max(size_v);

N_traces = length( sigchannels ) * length( trace_subset );
d0_signal = zeros( max_size, N_traces );
d0_reference_ladder = zeros( max_size, N_traces );


% If more than one color channel is specified, first fill a matrix with channel 1, then
%  fill additional traces from channel 2, in the same order...
count = 0;
tic
for m = 1:length(sigchannels) % usually just channel 1.
  for i = 1:length( trace_subset )  

    count = count+1;
    L = size(data_all{trace_subset(i)},1);
    % this is a straightforward subtraction of an offset.
    d0_signal          (1:L,count) = baseline_subtract(data_all{ trace_subset(i) }(:,sigchannels(m)));
    d0_reference_ladder(1:L,count) = baseline_subtract(data_all{ trace_subset(i) }(:,refchannel));

    labels{count} = filenames_all{ trace_subset(i) };
  end

  subset_pos = [subset_pos, count + subset_pos];
end
toc

% decide what labels to use
if isempty(lane_names); labels_used = labels; else; labels_used = lane_names; end;

% If user has not specified ymin,ymax in ylimit, figure it out.
AUTOFIND_YLIMIT = 0;
if isempty( ylimit ) 
  AUTOFIND_YLIMIT = 1;
  [ymin, ymax] = findTimeRange( d0_signal );
  ylimit = [ymin,ymax];
end
ymin = ylimit(1); ymax = ylimit(end);

if PLOT_STUFF
  h = figure(2); clf;
  set_print_page(h, 0, [50 50 600 800], 'All data');
  
  subplot(1,2,1);
  image( 0.05 * d0_signal);
  
  h=title( {dirnames{1}, sprintf('Signal channel(s) %s',  num2str( sigchannels) ) });
  set( h,'interpreter','none','fontweight','bold' );
  
  axis( [ 0.5 size( d0_signal, 2 )+0.5 ymin ymax] );
  set( gca, 'xtick', 1:size( d0_signal,2 ), ...
	    'xticklabel', char( labels_used )  );
  xticklabel_rotate;
  make_dividers( d0_signal, subset_pos, ymin, ymax );

  subplot(1,2,2);
  image( d0_reference_ladder*2);
  h=title( sprintf('Reference channel %s', num2str( refchannel) ));
  set( h,'interpreter','none','fontweight','bold' );   
  axis( [ 0.5 size( d0_signal, 2 )+0.5 ymin ymax] );
  make_dividers( d0_signal, subset_pos, ymin, ymax );
    set( gca, 'xtick', 1:size( d0_signal, 2), ...
	    'xticklabel', char( labels_used  )  );
  xticklabel_rotate;
  
  colormap( 1- gray(100));
end


% wipe out strong peaks from signal to leak into reference
if FIX_LEAKAGE_OF_SATURATING_SIGNALS_TO_REF
  for i = 1:length( data_all )
    data_all{i}(:,refchannel) = fix_leakage_of_saturating_signals_to_ref( data_all{i}(:,sigchannels(m) ),  data_all{i}(:,refchannel) );
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STAGE 3
reflane = 1; % this is a vestige of some old testing stuff.

group_count = 0;
for i = 1:length(data_set_starts)
  start_index = data_set_starts(i);
  if i < length( data_set_starts )
    final_index = data_set_starts(i+1)-1;
  else
    final_index = length( data_all );
  end
  if(start_index > final_index)
    continue;
  else
    group_count = group_count + 1;
  end
  data_align_group{group_count} = align_capillaries( ...
      { data_all{[start_index:final_index]} }, refchannel, reflane);

  ref_test = zeros(size(data_align_group{group_count}{1},1),length(data_align_group{group_count}));

  for j = 1:length(data_align_group{group_count})
    ref_test(:,j) = data_align_group{group_count}{j}(:,refchannel);
  end

  corr_test = zeros(size(ref_test,2));
  for j = 1:size(ref_test,2)
    for k = 1:size(ref_test,2)
      corr_test(j,k) = corr2(ref_test(:,j),ref_test(:,k));
    end
  end

  if(size(corr_test,2) > 2)
      [~,corr_idx] = sort(corr_test);

      corr_refined = zeros(size(corr_idx(2:end-1,:)));

      for j = 1:size(corr_refined,2)
          corr_refined(:,j) = corr_test(corr_idx(2:end-1,j),j);
      end

      corr_test = corr_refined;
  end

  corr_test = median(corr_test);

  bad_threshold = max(corr_test) / 2;

  if(corr_test(reflane) < bad_threshold)
    [~,another_ref] = max(corr_test);

    data_align_group{group_count} = align_capillaries( ...
          { data_all{[start_index:final_index]} }, refchannel, another_ref);
  end

end

data_align = align_capillaries_group( data_align_group, refchannel, 1, 1);


d = []; d_ref= [];

count = 0;
for m = 1:length(sigchannels) % usually just channel 1.

  for i = 1:length( trace_subset); 
    
    count = count+1;
    d(:, count)  = baseline_subtract(data_align{trace_subset(i)}(:, sigchannels(m)));
    d_ref(:, count) = abs(baseline_subtract(data_align{trace_subset(i)}(:, refchannel)));
    
    if(  ymax > size( d, 1 ) ) 
      fprintf( 'WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!\n')
      fprintf( 'WARNING! You specified a ymax (%5d) greater than trace length (%5d) WARNING!\n', ymax, size(d,1));
      ymax = size( d, 1);
      fprintf( 'WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!\n')
    end
    
    if NORMALIZE
      %d(:,i) = d(:,i)/mean( data_align{trace_subset(i)}(1500:2500,3) );
      d( :,count)    =  d(:,count)/mean( abs(d(ymin:ymax,count)));
      d_ref(:,count) = d_ref(:,count)/mean( abs(d_ref(ymin:ymax,count)));
    else
      d( :,count)    =  d(:,count);
      d_ref(:,count) = d_ref(:,count);
    end
    
  end
end;

if PLOT_STUFF
  h = figure(3); clf;
  set_print_page(h, 1, [100 100 600 800], 'Linear alignment');
  
  %subplot(1,2,1);
  image( 50*d );
  axis( [ 0.5 size( d, 2 )+0.5 ymin ymax] );
  set( gca, 'xtick', 1:size( d0_signal, 2 ), ...
	    'xticklabel', char( labels_used  )  );
  xticklabel_rotate;

  make_dividers( d, [], ymin, ymax );
  h=title( [dirnames{1}]);
  set( h,'interpreter','none','fontweight','bold' )
  colormap(  1 - gray(100) )

  %print( '-depsc2',[tag,'_Figure3.eps']);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STAGE 4 --  subtract off a smooth (but
%  not necessarily constant) baseline
figure(4);

ymin_pad = max(ymin - 200, 1);
ymax_pad = min(ymax + 200, size(d, 1));

if SMOOTH_BASELINE_SUBTRACT;   d = baseline_subtract_smooth( d, ymin_pad, ymax_pad);  end;

if (LOCAL_ALIGN) 
  ywindow = ymin_pad:ymax_pad;
  [ d(ywindow, :), d_ref(ywindow, :) ] = align_by_DP_using_ref( d(ywindow, :), d_ref(ywindow, :) ); 
end

d = d(ymin:ymax, :);
d_ref = d_ref(ymin:ymax, :);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STAGE 5
% align_by_DP -- local refinement through 
% piece-wise linear transformation.

if PLOT_STUFF
  h = figure(4); clf;
  set_print_page(h, 1, [150 150 600 800], 'Final (signal)');
 
  %subplot(1,2,1);
  image( 50*d );
  axis( [ 0.5 size( d0_signal, 2 )+0.5 1 size(d,1)] );
  %axis( [ 0.5 size( d0_signal, 2 )+0.5 ymin ymax] );
  set( gca, 'xtick', 1:size( d0_signal, 2 ), ...
	    'xticklabel', char( labels_used  ) ,'FontSize',8 );
  xticklabel_rotate;

  make_dividers( d0_signal, [], 1, size( d, 1 ));
  h=title( [dirnames{1}]);
  set( h,'interpreter','none' )
  colormap(  1 - gray(100) )

  %print( '-depsc2',[tag,'_Figure4.eps']);

  
  h = figure(5); clf;
  set_print_page(h, 1, [200 200 600 800], 'Final (reference)');
  
  %subplot(1,2,2);
  image( 50*d_ref );
  axis( [ 0.5 size( d0_signal, 2 )+0.5 1 size(d,1)] );
  %xlim( [ 0.5 size( d0_signal, 2 )+0.5 ] );
  set( gca, 'xtick', 1:size( d0_signal, 2 ), ...
	    'xticklabel', char( labels_used  )  );
  xticklabel_rotate;
  axis off
  title( 'aligned reference ladders');

  make_dividers( d0_signal, [], 1, size( d, 1 ));
  title( 'Reference ladder (channel 4)')
  
  colormap(  1 -gray(100) )
  
  figure(2)
  figure(4)
 
end

if PLOT_STUFF
  %fprintf( ['\nCreated: ',tag,'_Figure2.eps\n'] );
  %fprintf( ['Created: ',tag,'_Figure3.eps\n'] );
  %fprintf( ['Created: ',tag,'_Figure4.eps\n'] );
end


%% Clarence Cheng - save all figures as .eps and .fig
print_save_figure(figure(1), [tag,'_1Traces'], '', 1);
print_save_figure(figure(2), [tag,'_2AllData'], '', 1);
print_save_figure(figure(3), [tag,'_3LinearAlign'], '', 1);
print_save_figure(figure(4), [tag,'_4FinalSignal'], '', 1);
print_save_figure(figure(5), [tag,'_5FinalReference'], '', 1);


%%

figure(4);

fprintf( '\n' );
fprintf( 'ymin = %d\n', ymin)
fprintf( 'ymax = %d\n\n', ymax)
fprintf( '\n' );
fprintf( 'signal channel(s) = %s\n', num2str(sigchannels) )
fprintf( 'reference channel = %d\n\n', refchannel)
if ~isempty( dye_names_full ); 
  fprintf( 'Applied leakage correction for color channels:\n' ); 
  if ischar( dye_names_full );
    fprintf( [' ',dye_names_full,'\n'] );
  else
    for m = 1:length( dye_names_full ) fprintf( [' ',dye_names_full{m},'\n'] ); end;
  end
else  
  fprintf( 'No leakage correction applied to color channels.\n' ); 
end
if AUTOFIND_YLIMIT;              fprintf( 'Used auto-find of ymin, ymax.  [specify ylimit to turn off]\n' ); end;
if NORMALIZE;                    fprintf( 'Normalized data based on mean peak intensity. [set noNormalize to turn off]\n' ); end;
if SMOOTH_BASELINE_SUBTRACT;     fprintf( 'Applied subtraction of smooth base line [set noSmoothBaselineSubtract to turn off].\n' ); end;
if LOCAL_ALIGN;                  fprintf( 'Applied local alignment based on piece-wise linear transform [set noLocalAlign to turn off].\n' ); end;

fprintf( '\nFor all options, type: help %s\n', mfilename );

% gong notice when finished
gong

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  make_dividers( trace_subset, subset_pos, ymin, ymax )

hold on
for i = 1:size( trace_subset, 2 );
  plot( 0.5+i*[1 1], [ymin ymax],'k-', 'linew',0.25); 
end
for i = 2:length( subset_pos );
  plot( 0.5+subset_pos(i)*[1 1], [ymin ymax],'k-','linew',2); 
end
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dye_names_full = get_dye_names_full( dye_names, signals_and_ref )
dye_names_full = {};
if isempty( dye_names ); return; end;
if ischar( dye_names ); dye_names_full = dye_names; return; end;

%dye_names_full = {'FAM','HEX','TAMRA','ROX'};

for m = 1:length( signals_and_ref )
  dye_names_full{ signals_and_ref(m) } = dye_names{m};
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cols = split_string( l, delimiter )

if ~exist( 'delimiter','var') delimiter = ' '; end;

remain = l;
cols = {};
while ~isempty( remain )
  [token, remain] = strtok(remain, delimiter);
  cols = [cols, token];
end


%%
