function [d_out,x_transform_all, anchor_nodes] = align_by_DP_fine( d_all, align_blocks_in, dirnames )
% ALIGN_BY_DP_FINE: refine alignment by piece-wise-linear transform, optimizing correlation by dynamic programming
%
% calls ALIGN_BY_DP using parameter values that permit fine-grained realigment, as is useful
% with mutate/map data.
%
%  [d_out, x_transform_all, anchor_nodes] = align_by_DP_fine( d_all, align_blocks_in );
%
% Inputs:
%  d_all = matrix with traces to be aligned
%  align_blocks_in = subsets of traces for serial alignments specified as a cell of integer vectors. Example: 
%                       specifying { [1:4], [7 5 6 8] } will first align traces 2,3, and 4 to trace 1, and
%                         then trace 5, 6, and 8 to 7.  [Default is align all to 1]
%  
% Outputs:
%  d_out        = matrix with aligned traces
%  x_transform_all   = [advanced] matrix describing the local realignments
%  anchor_nodes = [advanced] window boundaries. 
%
% (C) R. Das, 2013
%

if nargin == 0;  help( mfilename ); return; end;

d_out = [];

if ~exist( 'align_blocks_in','var') align_blocks_in = []; end;

penalizeStretchFactor = 1.0;
slack = 10;
maxShift = 20;
windowSize = 50;
PLOT_STUFF = 1;
% maxShift changed from 10 to 20, by T47, tested on 16S mutate-and-map

[d_out, x_transform_all, anchor_nodes] = align_by_DP( d_all, align_blocks_in, penalizeStretchFactor, slack, maxShift, windowSize,  PLOT_STUFF );


% output to .eps file
if ~exist( 'dirnames','var') || isempty(dirnames);
    tag = '';
else
    tag = dirnames{1};
    if tag(end) == '/'; tag = tag(1:end-1); end; % get rid of final slash
    %  just get the last directory name.
    tags = split_string( tag, '/' );
    tag = tags{ end };
end;

if ~isempty(tag);
    if ~exist('Figures','dir'); mkdir('Figures'); end;
    tag = ['Figures/', tag, '_6DPAlign'];
    print( gcf, '-depsc2', '-loose', '-r300', [tag, '.eps']);
    fprintf( ['\nCreated: ', tag, '.eps\n'] );
    hgsave(gcf, tag);
end;
