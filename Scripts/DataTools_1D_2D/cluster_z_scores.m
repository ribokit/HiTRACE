function  cluster_z_scores( zscore_matrix, structure, offset  ); 
% cluster_z_scores( zscore_file, structure, offset  ); 
%
% (C) R. Das, 2011-2013.

if ischar( zscore_matrix ); 
  d = load( zscore_file );
else  
  d = zscore_matrix;
end

if ~exist( 'offset') offset = 0; end;
if ~exist( 'structure' ); structure = ''; end;

d = abs( smooth2d(d,1) );
figure(1)
image( d'*20 );
colormap( 1 - gray(100) );

set_axes( d, offset );

figure(2)
NRES = length( d ) ;
cluster_assignment = zeros( NRES, NRES );

count = 0;
CUTOFF = 1.0;
for i = 1:NRES
  for j = 1:NRES
    if d(i,j) > CUTOFF
      count = count + 1;
      cluster_assignment(i,j) = count;
    end
  end
end

cmap = jet( count );
cmap = cmap( randperm( count ), : );

cmap( 1,: ) = 0.0;
colormap( cmap );

NPASS = 5;
for n = 1:NPASS
  for i = 1:NRES
    for j = 1:NRES
      if ( cluster_assignment( i, j ) )
	cluster_assignment = check_neighbor( cluster_assignment, i, j, i-1, j );
	cluster_assignment = check_neighbor( cluster_assignment, i, j, i+1, j );
	cluster_assignment = check_neighbor( cluster_assignment, i, j, i,   j-1);
	cluster_assignment = check_neighbor( cluster_assignment, i, j, i,   j+1 );
	%cluster_assignment = check_neighbor( cluster_assignment, i, j, i-1, j+1 );
	%cluster_assignment = check_neighbor( cluster_assignment, i, j, i+1, j-1 );
	cluster_assignment = check_neighbor( cluster_assignment, i, j, j,   i );
      end
    end
  end
  

end

figure(1)
subplot(1,2,1)
image( cluster_assignment' );
set_axes( d, offset );
show_bps( structure );

% filter out weak clusters.
CLUSTER_CUTOFF = 8;
NUM_MUT_CUTOFF = 3; % disabled
ZSCORE_SUM_CUTOFF = 0.0; % disabled
MAX_ZSCORE_CUTOFF = 0.0; % disabled
MEAN_ZSCORE_CUTOFF = 0.0; % disabled
AVG_HITS_PER_MUT_CUTOFF = 1000; % disabled
CLUSTERS_PER_MUT_CUTOFF = 5; % important
FORCE_SYMM = 1;


cluster_assignment_plot = get_cluster_filter( d, cluster_assignment, CLUSTER_CUTOFF, ZSCORE_SUM_CUTOFF,MAX_ZSCORE_CUTOFF, MEAN_ZSCORE_CUTOFF, NUM_MUT_CUTOFF, AVG_HITS_PER_MUT_CUTOFF, FORCE_SYMM, CLUSTERS_PER_MUT_CUTOFF );

subplot(1,2,2)
image( cluster_assignment_plot' );
set_axes( d, offset );
show_bps( structure );


% for paper figure.
figure(2)
cmap( 1,: ) = 1.0;
colormap( cmap );
clf;
subplot(1,1,1);
image( cluster_assignment_plot' );
outline_clusters( cluster_assignment_plot, cmap );
set_axes( d, offset, 'k' );
show_bps( structure, 1 );
xticklabel_rotate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function show_bps( structure, USE_DOTS );

if length( structure ) == 0; return; end;
if ~exist( 'USE_DOTS' ); USE_DOTS = 0; end;

hold on;
bps = convert_structure_to_bps( structure );
bps = [ bps; bps(:,[2 1] ) ];
if USE_DOTS
  plot( bps(:,1), bps(:,2), 'ko','markersize',2,'markerfacecolor','w','linewidth',1 );
else
  for i = 1:size( bps, 1 )
    h = rectangle( 'Position', [ bps(i,1)-0.5, bps(i,2)-0.5, 1, 1], ...
		   'edgecolor', [0.7 0.7 0.7],...
		   'linewidth',0.5 );  
  end
end
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cluster_assignment  = get_cluster_filter( d, cluster_assignment, CLUSTER_CUTOFF, ZSCORE_SUM_CUTOFF, MAX_ZSCORE_CUTOFF, MEAN_ZSCORE_CUTOFF, NUM_MUT_CUTOFF, AVG_HITS_PER_MUT_CUTOFF, FORCE_SYMM, CLUSTERS_PER_MUT_CUTOFF );

NRES = length( d );

for i = 1:NRES
  n_clusters = length( unique( cluster_assignment( :, i ) ) );
  if ( n_clusters > CLUSTERS_PER_MUT_CUTOFF )
    cluster_assignment(:, i ) = 0;
  end
end

count = max( max( cluster_assignment) ) ;
for i = 1:count
  gp = find( cluster_assignment == i );
  if length(gp) > 0 & length(gp) < CLUSTER_CUTOFF
    cluster_assignment( gp ) = 0;
  end
end

for i = 1:count
  gp = find( cluster_assignment == i );
  if length(gp) > 0 & ( sum( d( gp ) ) < ZSCORE_SUM_CUTOFF | ...
			max( d( gp ) ) < MAX_ZSCORE_CUTOFF | ...
			mean( d(gp) ) < MEAN_ZSCORE_CUTOFF );
    cluster_assignment( gp ) = 0;
  end
end

[xgrid, ygrid ] = meshgrid( 1:NRES, 1:NRES );

for i = 1:count
  gp = find( cluster_assignment == i );
  if length(gp) > 0 
    num_mut_res = length( unique( xgrid( gp ) ) );
    if ( num_mut_res < NUM_MUT_CUTOFF )
      cluster_assignment( gp ) = 0;
    end
    if ( (length(gp) / num_mut_res ) >  AVG_HITS_PER_MUT_CUTOFF )
      cluster_assignment( gp ) = 0;
    end
    if FORCE_SYMM & ~check_for_symmetry( xgrid( gp ), ygrid( gp ) );
      cluster_assignment( gp ) = 0;
    end
  end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ok = check_for_symmetry( x, y )
ok = 0;
if ~isempty( intersect( x, y ) )
  ok = 1;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cluster_assignment = check_neighbor(  cluster_assignment, i, j, ix, jx );

if ix > size( cluster_assignment, 1); return; end;
if jx > size( cluster_assignment, 2); return; end;
if ix < 1 ; return; end;
if jx < 1 ; return; end;

if cluster_assignment( ix, jx )  
  cluster_assignment = replace_cluster_assignment( cluster_assignment, ...
						   cluster_assignment(i,j), ...
						   cluster_assignment(ix, jx ) );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function   cluster_assignment = replace_cluster_assignment( cluster_assignment, a, b );

if (a == b); return; end;
gp = find( cluster_assignment == a );
cluster_assignment( gp ) = b;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function set_axes( d, offset, axiscolor )
if ~exist( 'axiscolor' );  axiscolor = 'blue'; end;
nres = length( d );
seqpos = [1:nres] + offset;
gp = find( mod(seqpos,10) == 0 );
set( gca,'xtick', gp, 'xticklabel',  seqpos( gp ), 'ytick', gp, 'yticklabel', seqpos( gp ) );
set( gca,'xgrid','on','ygrid','on','fontsize',12,'fontweight','bold','xcolor',axiscolor,'ycolor',axiscolor );
hold on
plot( [1:nres],[1:nres],axiscolor);
hold off
set(gcf,'color','white');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outline_clusters( cluster_assignment_plot, cmap );

hold on

nres = size( cluster_assignment_plot, 1);
for i = 1:nres
  for j = 1:nres
    c = cluster_assignment_plot(i,j);
    if ( c > 0 )
      if ( i == 1 | ~cluster_assignment_plot(i-1,j)  )
	plot( [i-0.5 i-0.5], [j-0.5 j+0.5], 'color', 0.2 * cmap( c,: ) );
      end
      if ( i ==nres  | ~cluster_assignment_plot(i+1,j)  )
	plot( [i+0.5 i+0.5], [j-0.5 j+0.5], 'color', 0.2 * cmap( c,: ) );
      end
      if ( j == 1 | ~cluster_assignment_plot(i,j-1)  )
	plot( [i-0.5 i+0.5], [j-0.5 j-0.5], 'color', 0.2 * cmap( c,: ) );
      end
      if ( j == nres | ~cluster_assignment_plot(i,j+1)  )
	plot( [i-0.5 i+0.5], [j+0.5 j+0.5], 'color', 0.2 * cmap( c,: ) );
      end
    end
  end
end

hold off