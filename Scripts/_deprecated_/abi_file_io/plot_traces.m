function [a] = plot_traces( sarray, refarray, lanes, xrange, legends, plottitle, refplot, NORM, newfig, filename)

%   Use this function to plot lane traces from ABI sequencer data and analyzed by HiTRACE.
%   Read data into an array using the quick_look function.
%       (e.g. [x,xref]=quick_look(...), or [x,xref]=quick_look2(...), or [x,xcorr,xlabels]=quick_lookMC(...) to save quick_look'd data)

%Inputs:
%   sarray:     Name of array containing analyzed data (e.g. x)**
%   refarray:   Name of array containing analyzed reference data (e.g. xref)**
%   lanes:      Specify the lanes for which to view traces (e.g. [1,2:4])
%   xrange:     Specify the minimum and maximum time points to plot (e.g. [300 3500])       % need to edit so that entering a range (e.g. 300:3500) will prompt script to print those values as x-axis
%   legends:    Enter a list of legend names (e.g. {'leg1','leg2','leg3',...}, otherwise enter 'off'
%   plottitle:  Enter a string title
%   refplot (option):   Enter 1 to plot reference trace; 0, '', or omit otherwise
%   normalize (option): Enter 1 to normalize signals across lanes; 0, '', or omit otherwise
%   newfig (option):    Enter 1 to open in new figure; 0, '', or omit otherwise
%   filename (option):  Saves plot as PostScript and MatLab figure under given string filename, omit for no save
%
%   **if calling a cell array, as those produced by quick_look with multiple color channels, reference the cell corresponding to the desired
%     channel (e.g. x{1,1} for channel 1, x{1,2} for channel 2, etc)
%
%   Examples: plot_traces(xarray,xarray_ref,1:5,[1 3000],{'label1','label2'},'Plot Title',0,0,0,'Filename');
%             plot_traces(xarray,xarray_ref,:,[500 5000],'off','Plot Title');
%
%(C) Clarence Cheng, 2012

if ~exist( 'refplot' ) || length(refplot) == 0; refplot = 0;end;
if ~exist( 'newfig' ) || length(newfig) == 0; newfig = 0;end;
if ~exist( 'NORM' ) || length(NORM) == 0; NORM = 0;end;

if newfig == 1
    figure;
elseif newfig == 0
end

%Optionally normalize each lane to its max signal
if NORM == 1
    for i = [1:length(lanes)]
        sig_max = max(sarray(:,lanes(i)));
        normalize(:,i) = sarray(:,lanes(i))/sig_max * 100;   %build new array with normalized signals
    if refplot == 1         %Optionally plot reference trace with each sample trace
        hold on
        plot(normalize(:,i));
        plot(refarray(:,lanes),'--r');
        ylim([0 110]);
        hold off
    elseif refplot == 0
        plot(normalize);
    end
    end
elseif NORM == 0
    if refplot == 1         %Optionally plot reference trace with each sample trace
        hold on
        plot(sarray(:,lanes));          %xrange(1):xrange(2),
        plot(refarray(:,lanes),'--r');
        hold off
    elseif refplot == 0
        plot(sarray(:,lanes))
    end    
end

%Set range of x (time) values
xlim(xrange);


%Display legend with desired text in upper-LH corner
legend(legends, 'Location', 'SouthOutside','Orientation','horizontal');

%Display desired title
title(plottitle);


%Label x and y axes
xlabel('Time (ms)');
ylabel('Intensity');


%Save as postscript and figure under filename
if ~exist('filename')
else
    print(filename,'-dpsc');
    hgsave(filename);
end