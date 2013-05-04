function perm_mod = make_modifier_dendrogram( reactivity_final, subset_seq, subset_mod, seqpos, tags_final );
% MAKE_MODIFIER_DENDROGRAM
%
% perm_mod = make_modifier_dendrogram( reactivity_final, subset_seq, subset_mod, seqpos, tags_final );
%
%  Clustering of deep chemical profiles and a nice big plot.
%
clf;

%r = quick_norm( max(reactivity_final( subset_seq,subset_mod),0) );
r = quick_norm( max(remove_offset(reactivity_final( subset_seq,subset_mod) ),0) );

dm = pdist( r' ,'correlation' );
%dm = pdist( r' ,'mydist' );

z = linkage( dm, 'weighted' );
[h,t,perm_mod] = dendrogram( z, 0, 'labels',tags_final( subset_mod) );
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
