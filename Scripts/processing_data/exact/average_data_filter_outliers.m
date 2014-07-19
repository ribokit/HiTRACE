function [final_reactivity, final_error, flags ] = average_data_filter_outliers( reactivities, guessed_errors, seqpos, sequence, offset, trace_legends ); 
% AVERAGE_DATA_FILTER_OUTLIERS
% [final_reactivity, final_error, flags ] = average_data_filter_outliers( reactivities, guessed_errors ); 
%
% Takes traces and initial error estimates; figures out outlier points and even outlier traces;
%  and then returns reasonable final values and error estimates.
%
% Inputs:
%  reactivities   = Matrix of input traces
%  guessed_errors = Matrix of corresponding error estimates. These are used to weight the average and
%                     filter outliers, but final_error depends on the actual scatter seen in (non-outlier) points.
%  seqpos         = [Optional] conventional sequence numbering -- for warnings.
%  sequence       = [Optional] sequence -- for plotting
%  offset         = [Optional] offset to get from 1:N to conventional numbering -- for plotting
%  trace_legends  = [Optional] cell of legends for each trace
%
%
% Outputs:
%  final_reactivity = final averaged reactivity, weighted by 1/error^2.
%  final_error      = error on final reactivity, from standard deviation across traces / sqrt(N).
%  flags            = matrix of 0 and 1's, with 1 flagging problem points.
%
% (C) R. Das, T. Mann, Stanford University, 2013.

POINT_ERROR_RATIO_CUTOFF = 5.0;
PROFILE_ERROR_RATIO_CUTOFF = 2.5;

if ( nargin < 2 ) help( mfilename ); return; end;

N = size( reactivities, 2 );
L = size( reactivities, 1 );
if ~exist( 'seqpos','var') seqpos = [1:L]; end;
if ~exist( 'sequence','var') sequence = ''; end;
if ~exist( 'offset', 'var' ) offset = 0; end;
if ~exist( 'trace_legends', 'var' ) trace_legends = {}; end;

flags = ones(L,N);
final_err = zeros(L,1);
final_reactivity = zeros(L,1);

NITER = 5;
for i = 1:L
  gp = 1:N;

  weights = max( 1./guessed_errors(i,:), 0 ) ;
  for n = 1:NITER
    m = sum( reactivities( i, gp ) .*weights(gp) ) / sum( weights(gp) ) ;
    dev = abs(reactivities(i,:) - m );
    gp = find( dev < POINT_ERROR_RATIO_CUTOFF * guessed_errors(i,:)  );
  end  
  flags(i,gp) = 0;
  
  if length( gp ) == 0; % whoa that's a big problem.
    fprintf( 'Residue %d has way more scatter than input guessed_errors!\n', seqpos(i) );
    gp = [1:N];
  end
  
  m = sum( reactivities( i, gp ) .* weights(gp) ) / sum( weights(gp) );
  final_reactivity(i) = m; 
  
  if length( gp ) < 2; gp = [1:N]; end;
  s = std( reactivities(i, gp) );
  final_error(i) = max( s/sqrt( length(gp)), 0 );
end
final_error = final_error';

% now look over all input reactivity profiles -- are there some that are just crazy?
bad_traces = [];
for j = 1:N
  mean_rel_error(j) = mean( max(abs(reactivities(:,j) - final_reactivity) ./ final_error,0) );
end
for j = 1:N  
  %[mean_rel_error(j) median( mean_rel_error )]
  if ( mean_rel_error(j) > (PROFILE_ERROR_RATIO_CUTOFF * median( mean_rel_error )) )
    bad_traces = [bad_traces, j];
  end
end
good_traces = setdiff( [1:N], bad_traces );

if length( good_traces ) > 3 & length( bad_traces ) > 0
  fprintf( ['\nFollowing traces look way off -- will not use them for averaging: ',num2str(bad_traces),'\n\n'] );
  flags(:,bad_traces ) = 1.0;

  good_traces = setdiff( [1:N], bad_traces );
  [final_reactivity, final_error, flags(:,good_traces)] = average_data_filter_outliers( reactivities(:,good_traces), guessed_errors(:,good_traces), seqpos );
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

if length( trace_legends ) > 2; 
  h = legend( trace_legends, 2 ); 
  set(h,'interp','none','fontsize',6,'position',[0.3 0.4 0.05 0.05]);
end
draw_sequence( seqpos, sequence, offset, -0.1 );
ylabel('Reactivity');
ylim( [-0.5 4] );

set(gcf, 'PaperPositionMode','auto','color','white');
drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function draw_sequence( seqpos, sequence, offset, seq_level );
set(gca,'fontsize',12,'fontweight','bold','tickdir','out');
if length( sequence ) > 0
  for i = 1:length( seqpos ); 
    text( seqpos(i), seq_level,sequence(seqpos(i)-offset),'verticalalign','top','fontsize',7,'horizontalalign','center' ); 
  end
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
  plot( seqpos(i)*[1 1], r(i) + e(i)*[-1 1], 'color',c,'linewidth',lw );
  plot( seqpos(i) + [-0.5 0.5], r(i) + e(i)*[1 1], 'color',c,'linewidth',lw );
  plot( seqpos(i) + [-0.5 0.5], r(i) - e(i)*[1 1], 'color',c,'linewidth',lw );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function draw_flags_as_squares( flags, seqpos );
L = size( flags, 1 );
N = size( flags, 2 );
for i = 1:L
  for j = 1:N
    if ( flags(i,j) )
      rectangle( 'position', [seqpos(i)-0.5, j-0.5, 1, 1], 'edgecolor','r','linew',1.5 );
    end    
  end
  if all(flags(i,:) )
    rectangle( 'position', [seqpos(i)-0.5, -0.5, 1, N+1 ], 'edgecolor','r','linew',2 );
  end
end



