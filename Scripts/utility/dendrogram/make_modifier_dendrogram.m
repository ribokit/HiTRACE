function perm_mod = make_modifier_dendrogram( reactivity_final, subset_seq, subset_mod, seqpos, tags_final );

clf;

r = quick_norm( max(reactivity_final( subset_seq,subset_mod),0));
dm = pdist( r' ,'correlation' );
z = linkage( dm, 'weighted' );
[h,t,perm_mod] = dendrogram( z, 0, 'labels',tags_final( subset_mod) );
set(gca,'Position', [0.13 0.65 0.775 0.33] )
xlim([0.5 length(perm_mod)+0.5]);
xticklabel_rotate

subplot(2,1,2);
image([1:length(perm_mod)], seqpos(subset_seq), 40*r(:,perm_mod) )
box off
set(gca,'tickdir','out');
make_lines([1:length(perm_mod)] );
colormap( 1 - gray(100) );
set(gcf, 'PaperPositionMode','auto','color','white');
%print( '-depsc2','RhijuWinners_dendrogramALL_modifiers.eps');
