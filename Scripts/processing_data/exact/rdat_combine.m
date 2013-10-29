function [ final_rdat,flags ] = rdat_combine( rdat_files, outfilename)
% RDAT_COMBINE - 
% [ final_rdat,flags ] = rdat_combine( rdat_files, outfilename)
%
% Combines multiple rdat files into one, whose reactivity
%   values come from an average of all profiles, weighted by uncertainty.
% Final error values come from standard deviation across traces [divided by sqrt(N)],
%   with outliers ignored
% Flags positions or entire traces that agree poorly between replicates.  Assumes that
%  all the data is the same between input rdats other than the filename,
%  reactivity, and reactivity_error.
%
%
% Inputs:
%   rdat_files  = cell of strings containing rdat filenames to combine
%   outfilename = filename for final RDAT, formatted for RMDB
%
% Outputs:
%   final_rdat  = unified .rdat file with weighted average reactivities
%   flags       = residue positions that agree poorly with others
%
% (C) T. Mann, R. Das, Stanford University, 2013.
%

if ( nargin < 1 | length( rdat_files ) < 1 ); help( mfilename); return; end;
rdats = {};
reactivities = {}; errors = {}; final_reactivity = {}; prop_err = {};

%reads in all reactivities and errors together
N_rdat = length( rdat_files );
MIN_REL_ERROR = 0.1;
N = 0;
for i = 1:N_rdat
  rdats{i} = read_rdat_file(rdat_files{i});

  reactivity = rdats{i}.reactivity;
  for j = 1:size( reactivity,2 )
    if all( reactivity(:,j) == 0 ); continue; end;
    N = N + 1;
    reactivities{N} = reactivity(:,j);
    min_error = mean( max(reactivity(:,j),0) ) * MIN_REL_ERROR;
    errors{N} = max( rdats{i}.reactivity_error(:,j), min_error );
    rdat_file_legends{N} = basename( rdat_files{i} );
  end
  
end;

for i = 1:N
  for j = 1:size( reactivities{i} )
    if ( reactivities{i}(j) < 0 ) errors{i}(j) = max( errors{i}(j), abs( reactivities{i}(j) ) ); end;
  end
end

reactivities = cell2mat(reactivities);
errors = cell2mat(errors); %reformatting to use std function
L = size( reactivities, 1 );
sequence = rdats{1}.sequence;
offset = rdats{1}.offset;
seqpos = rdats{1}.seqpos;

[final_reactivity, final_error, flags ] = average_data_filter_outliers( reactivities, errors, seqpos ); %averages data, weighted by uncertainty

final_error = final_error;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make plots for visual feedback.
set(gcf,'position',[ 10  300  1200  500]);
clf;
subplot(2,1,1);
colormap( 1 - gray(100))
set(gca,'Position',[0.05 0.6 0.92 0.35] );
image( seqpos, [1:N], 128 * reactivities' );
make_lines_horizontal( [0:N],'k',0.5 );
draw_sequence( seqpos, sequence, offset, N+0.5 );
ylabel('Replicates');
draw_flags_as_squares( flags, seqpos );
h=title( outfilename );set(h,'interp','none','fontweight','bold');
xlabel('Sequence Position');

subplot(2,1,2);
set(gca,'Position',[0.05 0.1 0.92 0.35] );
colorcode = [ 0 0 1; 0 0.5 0; 0.5 0 0; 1 0 1; 1 0.5 0; 0.5 0.5 0.5; 0.5 0.5 1; 0.2 0.7 0.2; 1 0.5 0.5; 0.5 0 0.5; 1 0.7 0.3 ];
if size( colorcode, 1 ) < N; colorcode = [ colorcode; jet( N-size(colorcode,1) ) ]; end;
for j = 1:N;  make_plot( seqpos, reactivities(:,j), colorcode(j,:), 1 ); hold on; end;
plot( seqpos, 0*seqpos+1,'-','color',[0.7 0.7 0.7] );
make_plot(seqpos,final_reactivity,'k',1.5);
for j = 1:N; make_flags( seqpos, reactivities(:,j), flags(:,j) );end;
%  for i = 1:N;  make_plot_errors( seqpos, reactivities(:,i), errors(:,i), colorcode(i,:) ); hold on; end;
make_plot_errors(seqpos,final_reactivity,final_error,'k',1);
plot( seqpos, 0*seqpos, 'k' );

h = legend( rdat_file_legends, 2 ); 
set(h,'interp','none','fontsize',6,'position',[0.3 0.4 0.05 0.05]);
draw_sequence( seqpos, sequence, offset, -0.1 );
ylabel('Reactivity');
ylim( [-0.5 4] );

num_flags = [];
set(gcf, 'PaperPositionMode','auto','color','white');
drawnow;

if ~exist( 'Figures', 'dir' ) mkdir( 'Figures/' ); end;
export_fig(  ['Figures/',basename(outfilename), '.pdf'] );
save( ['Figures/',basename(outfilename), '.fig'] );

%RDAT PREPARATION using new filename, reactivity, and reactivity_error
name = rdats{1}.name;
if isempty(rdats{1}.structure); structure = []; else structure = rdats{1}.structure; end;
if isempty(rdats{1}.annotations); annotations = {}; else annotations = rdats{1}.annotations; end;
if isempty(rdats{1}.data_annotations); data_annotations = {}; else data_annotations = rdats{1}.data_annotations; end;
if isempty(rdats{1}.trace); trace_in = []; else trace_in = rdats{1}.trace_in; end;
if isempty(rdats{1}.xsel); xsel = []; else xsel = rdats{1}.xsel; end;
if isempty(rdats{1}.xsel_refine); xsel_refine = []; else xsel_refine = rdats{1}.xsel_refine; end;
if isempty(rdats{1}.comments); comments = []; else comments = rdats{1}.comments; end;

final_rdat = output_workspace_to_rdat_file( outfilename, name, sequence, offset, ...
			       seqpos, final_reactivity , ...
			       structure, ...
			       annotations, data_annotations, final_error,...
			       trace_in,xsel,xsel_refine,comments );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function draw_sequence( seqpos, sequence, offset, seq_level );
set(gca,'fontsize',12,'fontweight','bold','tickdir','out');
for i = 1:length( seqpos ); 
  text( seqpos(i), seq_level,sequence(seqpos(i)-offset),'verticalalign','top','fontsize',7,'horizontalalign','center' ); 
end
box off
xlim( [min(seqpos)-1 max( seqpos )+1] );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  make_plot( seqpos, r, c, lw );
if ~exist( 'lw','var') lw = 1; end;
if ~exist( 'c','var') c = 'k'; end;
plot( seqpos, r, '-','markerfacecolor',c,'color', c, 'linewidth',lw );
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  make_flags( seqpos, r, flags, c, lw );
if ~exist( 'lw','var') lw = 1; end;
if ~exist( 'c','var') c = 'r'; end;
hold on
if all( flags );  % distracting show circles if the whole thing is an outlier trace.
  plot( seqpos, r, 'x','color',c,'markersize',5 );
  return
end; 
bp = find( flags );
plot( seqpos(bp), r(bp), 'o','color',c );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  make_plot_errors( seqpos, r, e, c, lw );
if ~exist( 'lw','var') lw = 1; end;
if ~exist( 'c','var') c = 'k'; end;
hold on
for i = 1:length(seqpos )
  plot( seqpos(i)*[1 1] + 0.02*randn(1), r(i) + e(i)*[-1 1], 'color',c,'linewidth',lw );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function draw_flags_as_squares( flags, seqpos );
L = size( flags, 1 );
N = size( flags, 2 );
for i = 1:L
  for j = 1:N
    if ( flags(i,j) )
      rectangle( 'position', [seqpos(i)-0.5, j-0.5, 1, 1], 'edgecolor','r' );
    end    
  end
  if all(flags(i,:) )
    rectangle( 'position', [seqpos(i)-0.5, -0.5, 1, N+1 ], 'edgecolor','r','linew',2 );
  end
end



%%%%%%%%%%%%%%%
function b = basename( tag );

remain = tag;

while ~isempty( remain )
  [token, remain ] =strtok( remain, '/' );
end

b = token;