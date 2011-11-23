function [d,da,labels,d_noalign,da_noalign] = quick_look( dirnames, ymin, ymax, reorder, ...
				     labels, refcol, PLOT_STUFF, reflane )
% QUICK_LOOK  basic read-in script for ABI sequencer files for HiTRACE analysis
%
%  [d,da,labels,d_noalign,da_noalign] = quick_look( dirnames, ymin, ymax, reorder, labels, refcol, PLOT_STUFF )
%
%   quick_look() runs scripts to read in directories of .ab1 files from ABI 
%      sequencers and align them via a reference channel, assumed by default to be 
%       channel 4 of 4 (the texas red channel in Das lab experiments.) The returned output 
%      is d (aligned signal channel -- assumed to be channel 1) and da (aligned reference channel).
%  
%  Required input is:
%
%     dirnames = list of directories names  {'dir1','dir2',...'dirN'}
%
%  The rest of the inputs are optional.
%
%     ymin = minimum time point to show (e.g., 500, only used for plotting, script returns full profiles)
%     ymax = maximum time point to show (e.g., 3500, only used for plotting, script returns full profiles)
%     reorder = subset of electropherograms to align and return (e.g., [1:8 12:16])
%     labels = any user-defined labels (if you give as [], script will use filenames as labels on plots)
%     refcol = reference channel (Default 4. Other common setting would be the green channel 1 )
%     PLOT_STUFF = setting used by GUI interface to turn off plots.
%
%
% (C) Rhiju Das, 2009-2011
%
if ~exist( 'ymin' ); ymin = 500;   ymax = 3500; end
if exist( 'PLOT_STUFF' ) ~= 1;  PLOT_STUFF = 1; end;
if ~exist( 'reflane' ); reflane = 1;end;

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

line_pos = [ 0 ];
filepath = '';
[ data_all, filenames_all, data_init, data_length ] = ...
    read_abi_dirs( filepath, dirnames, PLOT_STUFF );

if length( data_all ) == 0 ; return; end;

line_pos = data_init-1;
numfiles = length( data_all );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STAGE 2
%Quick check... let's plot the raw data.
% The signal we care about is in channel 1.
% A reference ladder (in a different color) shows up in channel 4.
if ~exist( 'reorder' ) | length( reorder) == 0;  reorder = [ 1 : length( data_all ) ]; end

% Copy the appropriate columns to a single matrix, to make plotting easy.
% standard...
if ~exist( 'refcol');  refcol = [ 4 ]; end

%truncate data to the same length (this shouldn't be necessary, but happened once with a mixed-up data set from the PAN facility.
%for i = 1:length( data_all )
%  if (i == 1); whichpixels = [ 1 : size( data_all{reorder(i)}, 1 ) ]; end;
%  data_all{i} = data_all{i}(whichpixels,:);
%end

for i = length(reorder)
    size_v = size(data_all{i},1);
end
max_size = max(size_v);

d0_signal = [];
d0_reference_ladder = [];

%for i = 1:length( reorder )  
%  d0_signal(:,i)           = baseline_subtract(data_all{reorder(i)}(:,1));
%  d0_reference_ladder(:,i) = baseline_subtract(data_all{reorder(i)}(:,refcol));
%end
for i = 1:length( reorder )  
    if(size(data_all{reorder(i)},1) < max_size)
        data_all{reorder(i)}(max_size, size(data_all{reorder(i)},2)) = 0;
    end
    
  d0_signal(1:size(data_all{reorder(i)}(:,1),1),i)           = baseline_subtract(data_all{reorder(i)}(:,1));
  d0_reference_ladder(1:size(data_all{reorder(i)}(:,1),1),i) = baseline_subtract(data_all{reorder(i)}(:,refcol));
end

%for output, in case we want it.
d_noalign = d0_signal;
da_noalign = d0_reference_ladder;

if ~exist('labels') | length( labels ) == 0 
  labels = filenames_all;
end

if PLOT_STUFF
  figure(2)
  clf
  subplot(1,2,1);
  image( 0.05 * d0_signal);
  
  h=title( [dirnames{1}]);
  set( h,'interpreter','none' );
  
  axis( [ 0.5 length( reorder )+0.5 ymin ymax] );
  set( gca, 'xtick', 1:length( reorder ), ...
	    'xticklabel', char( labels{reorder}  )  );
  xticklabel_rotate;
  
  for i = 1:length( line_pos );
    hold on
    plot( 0.5+line_pos(i)*[1 1], [ymin ymax],'r-'); 
  end
  hold off


  %figure(3)
  subplot(1,2,2);
  image( d0_reference_ladder*2);
  title( 'Reference ladder (channel 4)')
  axis( [ 0.5 length( reorder )+0.5 ymin ymax] );
  
  for i = 1:length( line_pos );
    hold on
    plot( 0.5+line_pos(i)*[1 1], [ymin ymax],'r-'); 
  end
  hold off
  
  set( gca, 'xtick', 1:length(reorder), ...
	    'xticklabel', char( labels{reorder}  )  );
  xticklabel_rotate;
  
  colormap( 1- gray(100));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STAGE 3

% refine?
%data_align = align_capillaries( data_all, [4], 1, 1 );
%data_align = align_capillaries( data_align, [4], 1 );

%data_align = align_capillaries( data_all, refcol, 1 );
%data_align = align_capillaries( data_align, [1], 1, 1 );

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


NORMALIZE = 1;
for i = 1:length( reorder); 
  d(:,i)  = baseline_subtract(data_align{reorder(i)}(:,1));
  da(:,i) = baseline_subtract(data_align{reorder(i)}(:,refcol));

  if(  ymax > size( d, 1 ) ) 
    fprintf( 'WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!\n')
    fprintf( 'WARNING! You specified a ymax (%5d) greater than trace length (%5d) WARNING!\n', ymax, size(d,1));
    ymax = size( d, 1);
    fprintf( 'WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!\n')
  end

  if NORMALIZE
    %d(:,i) = d(:,i)/mean( data_align{reorder(i)}(1500:2500,3) );
    d(:,i) =  d(:,i)/mean( abs(d(ymin:ymax,i)));
    da(:,i) = da(:,i)/mean( abs(da(ymin:ymax,i)));
  else
    d(:,i) = d(:,i)/100;
    da(:,i) = da(:,i)/50;
  end
  
end;

%save d.mat d;
%save da.mat da;

if PLOT_STUFF
  figure(3)
  clf
  
  %subplot(1,2,1);
  image( 50*d );
  axis( [ 0.5 length( reorder )+0.5 ymin ymax] );
  set( gca, 'xtick', 1:length( reorder ), ...
	    'xticklabel', char( labels{ reorder}  )  );
  xticklabel_rotate;
  
  hold on
  for i = 1:length( line_pos );
    hold on
    %plot( 0.5+line_pos(i)*[1 1], [ymin ymax],'r-'); 
  end
  hold off
  h=title( [dirnames{1}]);
  set( h,'interpreter','none' )
  colormap(  1 - gray(100) )
  
  figure(4)
  %subplot(1,2,2);
  image( 50*da );
  axis( [ 0.5 length( reorder )+0.5 ymin ymax] );
  set( gca, 'xtick', 1:length( reorder ), ...
	    'xticklabel', char( labels{ reorder }  )  );
  xticklabel_rotate;
  axis off
  title( 'aligned reference ladders');
  
  hold on
  for i = 1:length( line_pos );
    hold on
    plot( 0.5+line_pos(i)*[1 1], [ymin ymax],'r-'); 
  end
  hold off
  title( 'Reference ladder (channel 4)')
  
  colormap(  1 -gray(100) )
  
  figure(3)
  
  labels = labels( reorder );
end