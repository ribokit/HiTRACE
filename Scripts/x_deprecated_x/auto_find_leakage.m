function [lm,d,data_bsub] = auto_find_leakage(dirnames,ymin,ymax,reorder, xmin,xmax, inlm, PLOT_STUFF);
% Function to give the leakage matrix between four channels returned from
% ABI sequencer. First load subsaturating amounts of each of four dyes into
% four different lanes on ABI capillary sequencer. 
%
% lm = auto_find_leakage(dirnames,ymin,ymax,reorder, xmin,xmax, PLOT_STUFF);
%
% 'dirnames' is the directory of the four lanes each with one dye.
%
% 'ymin' and 'ymax' give the interval to display all four lanes.
%
% It is important to reorder the lanes such that the dye you want in
% channel one is first, the dye you want in channel two is second, etc. Use
% 'reorder' to set.
%
% The function goes through and integrates the signal from, for example, dye 1 in
% channels 2-4, normalized by the signal from dye 1 in
% channel 1. It can be advantageous to only integrate the signal around
% your band rather than the whole trace. Use 'xmin' and 'xmax' to set.
%
% 'xmin' and 'xmax' are vectors giving the desired integration intervals
% for Dye 1-4, i.e. xmin=[xmin_ch1 xmin_ch2 xmin_ch3 xmin_ch4]. By default,
% 'xmin' and 'xmax' are set to 'ymin' and 'ymax' for all four channels.



if ~exist( 'ymin' ); ymin = 1;   ymax = 5000; end
if exist( 'PLOT_STUFF' ) ~= 1;  PLOT_STUFF = 1; end;

%%


labels = {};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STAGE 1
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
%%

filepath = '';
[ data_all, filenames_all, data_init, data_length ] = ...
    read_abi_dirs( filepath, dirnames, PLOT_STUFF );

if length( data_all ) == 0 ; return; end;

numfiles = length( data_all );

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

%% Now determine leakage matrix.
% option to fit only from xmin to xmax. Otherwise, the fluoresence is fit
% from ymin to ymax.

if ~exist('xmin') || length(xmin)<4;
    xmin=ymin*ones(1,4);
end
if ~exist('xmax') || length(xmax)<4;
    xmax=ymax*ones(1,4);
end

lm=ones(4,4); % initialize leakage matrix

k=[1:4;5:8;9:12;13:16]; % used to call the lane in 16-column matrix d corresponding to each fluorophore/fluorescent channel

for i=1:4;
    lm(i,:)=sum(d(xmin:xmax,k(i,:)))/sum(d(xmin:xmax,k(i,i)));  %each fluorophore's fluorescence in one channel divided by the proper fluorophore's fluorescence in that channel
end
%%
d_correct=d;
for i=1:4;
    d_correct(:,k(i,:))=d(:,k(i,:))*lm^(-1);
end

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


%%
d_correct2=d;
for i=1:4;
    d_correct2(:,k(i,:))=d(:,k(i,:))*inlm^(-1);
end

%% Plot corrected gel
if PLOT_STUFF
  figure;
    figure(4)
  set(gcf, 'PaperPositionMode','auto','color','white');
  clf
  %subplot(1,2,1);
  image( d_correct2);
  
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

