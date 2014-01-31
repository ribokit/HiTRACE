function make_rdat_stair_plots( r, idx, colorcode, labels, norm_pos );
% make_rdat_stair_plots( r, idx, colorcode, labels, norm_pos );
%
% Inputs:
%
% r          = cell of rdats
% idx        = indices of reactivities within each rdat
% colorcode  = [Optional] RGB colors, e.g., [0 0 1; 1 0 0] will plot red, then blue. [Default: 'jet']
% labels     = [Optional] names for each rdat [default: 1, 2, 3, ... ]
% norm_pos   = [Optional] sequence positions over which to normalize.
%
% (C) R. Das, 2013

cla;
if nargin==0; help( mfilename ); return; end;

N = length( r );
set(gcf, 'PaperPositionMode','auto','color','white');

% plot values as staircase
if ~exist( 'idx', 'var' ) | isempty(idx); idx = ones( N, 1 ); end;
if ~exist( 'labels','var' ); for i = 1:N; labels{i} = sprintf( '(%d) %s', i, r{i}.name ); end; end;
if ~exist( 'colorcode' ) | isempty(colorcode); colorcode = jet( N ); end;
ymax = 0;

for i = 1:N; 

  if exist( 'norm_pos', 'var') & length( norm_pos ) > 0
    norm_bins = [];
    for m = norm_pos; norm_bins = [norm_bins, find( m == r{i}.seqpos ) ]; end;
    [r{i}.reactivity(:,idx(i)), dummy, r{i}.reactivity_error(:,idx(i))] = ...
		     quick_norm( r{i}.reactivity(:,idx(i)), norm_bins, r{i}.reactivity_error(:,idx(i)) );
  end

  
  stairs( r{i}.seqpos-0.5, ...
	  r{i}.reactivity( :, idx(i) ), ...
	  'color',colorcode(i,:), ...
	  'linew',2 ); 
  hold on; 
  ymax = max( ymax,  max( r{i}.reactivity( 4:end, idx(i) ) ) );
end

% plot error bars
for i = 1:N; 
  Nres = size( r{i}.reactivity, 1 );
  for j = 1:Nres; 
    val = r{i}.reactivity(j,idx(i));
    val_err = r{i}.reactivity_error(j, idx(i));
    plot( ( r{i}.seqpos(j))*[1 1], val + val_err*[1 -1], 'color',colorcode(i,:) );
  end
end

% make a line along x-axis
plot( [min( r{N}.seqpos)-0.5 max( r{N}.seqpos)+0.5],  [0 0], 'k' ); hold on; 
set(gca,'xaxisloc','top','xgrid','on','fontw','bold');

ymax = ymax*1.2;
if ymax>0; ylim([0 ymax]); end;

xmax = Nres;
last_sequence  =  r{N}.sequences{ idx(N) };
last_structure =  r{N}.structures{ idx(N) };
if length( last_sequence ) > 0; xmax = length( last_sequence ); end;
xlim( [-0.5 xmax+1] + r{N}.offset );
    
for j = 1:(Nres-1); 
  if ( j > length( last_sequence ) ); continue; end;
  seqchar = sprintf('%s', last_sequence(j) );
  if ( length( last_structure ) > 0 ) 
    seqchar = sprintf( '%s\n%s', seqchar,last_structure(j) );
  end
  h = text( j + r{N}.offset, 0, seqchar);
  set(h,'HorizontalAlignment','center','VerticalAlignment','top','fontw','bold');
end
hold off;

if ~exist( 'labels' ) labels = num2str( [1:N]' ); end;
for i = 1:length( labels ); labels{i} = strrep( labels{i}, '_', '\_' ); end;
h=legend( labels );