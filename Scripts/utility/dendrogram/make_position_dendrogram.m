function [perm_seq,perm_mod,dm_seq,dm_mod,r] = make_position_dendrogram( reactivity, reactivity_error, subset_seq, subset_mod, seqpos, labels, sequences, offsets, sources, source_names );
% MAKE_POSITION_DENDROGRAM
%
%  [perm_seq_out,perm_mod_out] = make_position_dendrogram( reactivity, subset_seq, subset_mod, seqpos, labels, sequences, offsets, sources, source_names );
%
% Still under testing.  
%
%

if ~exist( 'sources' ) sources = ones(1, length( seqpos ) ); end;
if ~exist( 'source_names' ) source_names = {}; end;
clf;
if ~iscell( sequences ); sequences = { sequences } ; end;
for i= 1:length(subset_mod); blank_tags{i} = ''; end;

% may want to put this back in...
%for i = 1:length( sources );
 % normpts = find( sources(subset_seq) == i );
 % r( normpts, :) = quick_norm( max(reactivity( subset_seq(normpts),subset_mod),0));
%end

%r = quick_norm( max(remove_offset(reactivity( subset_seq,subset_mod) ),0) );
r = max(remove_offset(reactivity( subset_seq,subset_mod) ), 0);
%r = 2*r/mean(mean(r));
[r, norm_factor, r_error] = quick_norm( r, [], reactivity_error( subset_seq, subset_mod) );


% new: cap_outliers
%r = min( max( r, 0 ), 5 );

%dm = pdist( r' ,'correlation' );
%dm = get_chi_squared_dm( r, r_error );

% this copies code from make_modifier_dendrogram. Hmm.
dm_mod = get_weighted_correlation_coefficient( r, r_error);
z = linkage( dm_mod, 'weighted' );
leaf_order_mod = optimalleaforder( z, dm_mod );
[h,t,perm_mod] = dendrogram( z, 0,'labels',blank_tags,'colorthreshold',0.3,'orientation','left','reorder',leaf_order_mod ); 

for k = 1:length(h); set(h(k),'linew',2); end;
set(gca,'Position', [0.01 0.05 0.05 0.80]);
ylim([0.5 length(perm_mod)+0.5] )
axis off
xticklabel_rotate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cluster over sequence positions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
seqpos_tags = {};
for i = 1:length(subset_seq); 
  m = seqpos(subset_seq(i)); 
  source = sources(subset_seq(i));
  seqpos_tags{i} = [sequences{ source }( m - offsets( source ) ),num2str(m) ]; 
  if length( source_names ) > 0; seqpos_tags{i} = [ source_names{source},':',seqpos_tags{i} ]; end;
end;


subplot(2,2,2);
%dm_seq = pdist( r ,'correlation' );
%dm_seq = pdist( r ,'mydist' );
%dm_seq = pdist( r ,'distcorr_wrapper' );
%dm_seq = get_chi_squared_dm( r', r_error' );
dm_seq = get_weighted_correlation_coefficient( r', r_error');

z_seq = linkage( dm_seq, 'weighted' );
%z_seq = linkage( dm_seq, 'ward' );
leaf_order_seq = optimalleaforder( z_seq, dm_seq );

[h,t,perm_seq] = dendrogram( z_seq, 0, 'labels',seqpos_tags,'colorthreshold',0.4,'reorder',leaf_order_seq );

for k = 1:length(h); set(h(k),'linew',2); end;
xlim([0.5 length(perm_seq)+0.5] )
set(gca,'Position', [0.15 0.9 0.80 0.1],'ticklength',[0 0],'tickdir','out','fontweight','bold','fontsize',9 );
xticklabel_rotate

subplot(2,2,4);
image(40 * r(perm_seq,perm_mod)' );
colormap( [1 - gray(100); zeros(200,3); 1 0 0] )
box off
set(gca,'tickdir','out' ,'ticklength',[0 0], 'ytick',[1:length(perm_mod)],'fontweight','bold','yticklabel',labels(subset_mod(perm_mod)),'fontsize',6 );
xlim( [0.5 length(perm_seq)+0.5] );
ylim( [0.5 length(perm_mod)+0.5] );
make_lines( [0:length(perm_seq)], 'k', 0.25);
make_lines_horizontal( [0:length(perm_mod)], 'k', 0.25);
set(gcf, 'PaperPositionMode','auto','color','white');
set(gca,'Position', [0.15 0.05 0.80 0.80],'ydir','normal')
set(gca,'xtick',1:length(perm_seq),'xticklabel',seqpos_tags(perm_seq),'fontweight','bold','fontsize',9 );
xticklabel_rotate;


perm_seq_out = subset_seq(perm_seq);
perm_mod_out = subset_mod(perm_mod);

set(gcf,'Pointer','fullcross');
