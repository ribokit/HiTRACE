function perm_mod = make_modifier_dendrogram( reactivity, reactivity_error, subset_seq, subset_mod, seqpos, labels );
% MAKE_MODIFIER_DENDROGRAM
%
% perm_mod = make_modifier_dendrogram( reactivity, subset_seq, subset_mod, seqpos, labels );
%
%  Clustering of deep chemical profiles and a nice big plot.
%
clf;

%r = quick_norm( max(reactivity( subset_seq,subset_mod),0) );

r = max(remove_offset(reactivity( subset_seq,subset_mod) ), 0);
[r, norm_factor, r_error] = quick_norm( r, [], reactivity_error( subset_seq, subset_mod) );

if ~isempty( find( sum( r )  == 0 ) )
  error( 'One of the input traces is zero!' );
end

% new: cap_outliers
%r = min( max( r, 0 ), 5 );

%dm = pdist( r' ,'correlation' );
%dm = pdist( r' ,'mydist' );
%dm = pdist( r' ,'distcorr_wrapper' );
dm = get_chi_squared_dm( r, r_error );

z = linkage( dm, 'weighted' );
[h,t,perm_mod] = dendrogram( z, 0, 'labels',labels( subset_mod) );
set(gca,'Position', [0.02 0.85 0.98 0.10] )
xlim([0.5 length(perm_mod)+0.5]);
axis off
xticklabel_rotate

subplot(2,1,2);
image([1:length(perm_mod)], seqpos(subset_seq), 40*r(:,perm_mod) )
box off
set(gca,'tickdir','in','ygrid','on');
make_lines([1:length(perm_mod)] );
colormap( 1 - gray(100) );
set(gcf, 'PaperPositionMode','auto','color','white');
set(gca,'Position', [0.02 0.02 0.98 0.63] )
%print( '-depsc2','RhijuWinners_dendrogramALL_modifiers.eps');

set(gcf,'Pointer','fullcross');
