function [lm, datacorr]=make_leakage_matrix(dirnames,ymin,ymax,reorder,peakmin,peakmax,PLOT_STUFF);
%
% Creates a leakage matrix to correct the signal detected by the four
% channels of an ABI3100 sequencer. Outputs this leakage matrix, as well as
% a leakage-corrected dataset.
%
% Load oligonucleotides tagged with the desired set of fluorophores in
% separate lanes; it is useful to include a series of dilutions to ensure
% that subsaturating peaks are observed.
%
% -'dirnames' is the directory of the four lanes each with one dye.
%
% -'ymin' and 'ymax' give the interval to display all four lanes.
%
% -It is important to reorder the lanes such that the dye you want in
%  channel one is first, the dye you want in channel two is second, etc. Use
%  'reorder' to set.
%
% -Peakmin and peakmax are user-defined ranges defining the locations of
%  the peaks for each fluorescently-labeled oligo used for the leakage
%  calibration experiment. Peakmin should be a vector with the minimum x
%  values of the peaks (ex. [min_ch1 min_ch2 min_ch3 min_ch4]) and peakmax
%  should contain the maximum x values of the peaks.
%
% -make_leakage_matrix integrates the user-defined peaks from each dye 
%
% Sarah Denny, 2012
% Clarence Cheng, 2013



if ~exist( 'ymin' ); ymin = 1;   ymax = 5000; end
if exist( 'PLOT_STUFF' ) ~= 1;  PLOT_STUFF = 1; end;

%%


labels = {};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STAGE 1
% The ABI file format is a special binary file format -- 
% Rhiju has written scripts to take care of the readin, and
% to put all the data into one object ("data_all")

if ~iscell( dirnames )
  if isdir( dirnames ) 
    dirnames = { dirnames };
  else 
    % Assume it is a text file with names.
    [dirnames] = textread(dirnames,'%s');
  end
end
%%

filepath = '';
[ data_all, filenames_all, data_init, data_length ] = ...
    read_abi_dirs( filepath, dirnames, PLOT_STUFF );

if length( data_all ) == 0 ; return; end;

numfiles = length( data_all );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STAGE 2
% Check that only four lanes are loaded and plot each lane in all four
% channels.
if ~exist( 'reorder' ) | length( reorder) == 0;  reorder = [ 1 : length( data_all ) ]; end

if length(reorder)<4; return; end;

for i = length(reorder)
    size_v = size(data_all{i},1);
end
max_size = max(size_v);

% Subtract linear baseline
data_bsub=data_all;

for i=1:length(reorder);
    if(size(data_all{reorder(i)},1) < max_size)
        data_all{reorder(i)}(max_size, size(data_all{reorder(i)},2)) = 0;
    end
    for j=1:4;
        data_bsub{i}(:,j) = baseline_subtract(data_all{reorder(i)}(:,j));
    end
end
% data_bsub is reordered 
% Put data_bsub into a matrix rather than a cell. Call it "d". For 4 lanes of data in 4 channels, it will be 16 columns long.
d = [];
for i=1:length(reorder)
    d=[d data_bsub{i}];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot uncorrected data
        
if ~exist('labels') | length( labels ) == 0 
  labels = filenames_all;

end
labels_real=cell(1,4*length(reorder));

for i=1:length(reorder)
    for j=1:4;
        labels_real{(i-1)*4+j}=[labels{reorder(i)} '_Ch' num2str(j)];
    end
end


if PLOT_STUFF
  figure(2)
  set(gcf, 'PaperPositionMode','auto','color','white');
  clf
  %subplot(1,2,1);
  image( d);
  
  h=title( 'Signal (channel 1)');
  set( h,'interpreter','none' );
  n=4*length(reorder);
  axis( [ 0.5 n+0.5 ymin ymax] );
  set( gca, 'xtick', 1:n, ...
	    'xticklabel', char( labels_real{1:n}  )  );
  xticklabel_rotate;
  line_pos=[4 8 12];
  for i = 1:length( line_pos );
    hold on
    plot( 0.5+line_pos(i)*[1 1], [ymin ymax],'r-'); 
  end
  hold off
  colormap( 1- gray(100));
  title('No leakage correction');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Now determine leakage matrix.
% option to fit only from peakmin to peakmax. Otherwise, the fluoresence is fit
% from ymin to ymax.

if ~exist('peakmin') || length(peakmin)<4;
    peakmin=ymin*ones(1,4);
end
if ~exist('peakmax') || length(peakmax)<4;
    peakmax=ymax*ones(1,4);
end

lm=ones(4,4); % initialize leakage matrix

k=[1:4;5:8;9:12;13:16];
% k is used to call the lane in 16-column matrix d corresponding to each
% fluorophore/fluorescent channel. Each row of k corresponds to all the
% lanes in one fluorescent channel, while each column of k corresponds to
% the signal of one fluorophore across all channels.

for i=1:4;
    lm(i,:)=sum(d(peakmin:peakmax,k(i,:)))/sum(d(peakmin:peakmax,k(i,i)));
    %each fluorophores' fluorescence in one channel divided by proper fluorophore's fluorescence in that channel
end
%%
d_correct=d;
for i=1:4;
    d_correct(:,k(i,:))=d(:,k(i,:))*lm^(-1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot corrected gel
if PLOT_STUFF
  figure(3)
  set(gcf, 'PaperPositionMode','auto','color','white');
  clf
  %subplot(1,2,1);
  image( d_correct);
  
  h=title( 'Signal (channel 1)');
  set( h,'interpreter','none' );
  n=4*length(reorder);
  axis( [ 0.5 n+0.5 ymin ymax] );
  set( gca, 'xtick', 1:n, ...
	    'xticklabel', char( labels_real{1:n}  )  );
  xticklabel_rotate;
  
  line_pos=[4 8 12];
  for i = 1:length( line_pos );
    hold on
    plot( 0.5+line_pos(i)*[1 1], [ymin ymax],'r-'); 
  end
  hold off
  colormap( 1- gray(200));
  title('Corrected for leakage');
end


