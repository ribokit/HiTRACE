function [perm_mod,dm,r] = make_modifier_dendrogram( reactivity, reactivity_error, subset_seq, subset_mod, seqpos, labels, seqpos_tags, USE_CORRELATION );
% MAKE_MODIFIER_DENDROGRAM
%
% [perm_mod,dm,r] = make_modifier_dendrogram( reactivity, reactivity_error, subset_seq, subset_mod, seqpos, labels, seqpos_tags );
%
%  Clustering of deep chemical profiles and a nice big plot.
%
% (C) R. Das, Stanford University, 2013.

if (nargin < 1); help( mfilename ); return; end;
  
if ~exist( 'seqpos_tags' ) | isempty( seqpos_tags ) seqpos_tags = [1:size( reactivity,1) ]; end;
if ~iscell( seqpos_tags ) & isnumeric( seqpos_tags ) seqpos = seqpos_tags; clear seqpos_tags; end;
if ~exist( 'USE_CORRELATION' ) USE_CORRELATION = 0; end;
clf;

[r, norm_factor, r_error ] = prepare_data( reactivity, subset_seq, subset_mod, reactivity_error );

if ~isempty( find( sum( r )  == 0 ) )
  error( 'One of the input traces is zero!' );
end

% new: cap_outliers
%r = min( max( r, 0 ), 5 );

%dm = pdist( r' ,'mydist' );
%dm = pdist( r' ,'distcorr_wrapper' );
%dm = get_chi_squared_dm( r, r_error );
if USE_CORRELATION
  dm = pdist( r' ,'correlation' );
else
  dm = get_weighted_correlation_coefficient( r, r_error );
end

%z = linkage( dm, 'ward' );
%z = linkage( dm, 'average' );
z = linkage( dm, 'weighted' );

leaf_order_mod = optimalleaforder( z, dm );
[h,t,perm_mod] = dendrogram( z, 0, 'labels',labels( subset_mod),'reorder',leaf_order_mod );
%[h,t,perm_mod] = dendrogram( z, 0,'reorder',leaf_order_mod );

set(gca,'Position', [0.05 0.85 0.95 0.10],'fontsize',7 )
xlim([0.5 length(perm_mod)+0.5]);
axis off
xticklabel_rotate

subplot(2,1,2);
image([1:length(perm_mod)], [1:length(subset_seq)], 40*r(:,perm_mod) )
box off
set(gca,'tickdir','in','ygrid','on','ticklength',[ 0 0] );
make_lines([1:length(perm_mod)] );
make_colormap;
set(gcf, 'PaperPositionMode','auto','color','white');
set(gca,'Position', [0.05 0.02 0.95 0.73] )
%print( '-depsc2','RhijuWinners_dendrogramALL_modifiers.eps');

set(gcf,'Pointer','fullcross');

dm = squareform(dm);

%if exist( 'seqpos_tags' ); 
set(gca,'ytick',1:length( subset_seq ), 'yticklabel', seqpos_tags( subset_seq ),'fontsize',6,'ygrid','off' );
gp = find( mod( seqpos, 10 ) == 0 );
make_lines_horizontal( seqpos(gp), 'k',0.25,':');
%else 
%  gp = find( mod(seqpos(subset_seq),10) == 0 );
%  set(gca,'ytick',subset_seq(gp),'yticklabel',seqpos(subset_seq(gp)) );
%end


