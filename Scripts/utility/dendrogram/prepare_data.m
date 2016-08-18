function [r, norm_factor, r_error ] = prepare_data( reactivity, subset_seq, subset_mod, reactivity_error );
% [r, norm_factor, r_error ] = prepare_data( reactivity, subset_seq, subset_mod, reactivity_error );
%
% Apply offset & normalization before making dendrograms or running t-SNE.
%
% (C) R. Das, 2016

if ~exist( 'reactivity_error', 'var' ); reactivity_error = 0 *reactivity; end;
if ~exist( 'subset_mod', 'var' ) || isempty( subset_mod ); subset_mod = [1:size(reactivity,2)]; end;

%r = quick_norm( max(reactivity( subset_seq,subset_mod),0) );
r = max(remove_offset(reactivity( subset_seq,subset_mod) ), 0);
[r, norm_factor, r_error] = quick_norm( r, [], reactivity_error( subset_seq, subset_mod) );
