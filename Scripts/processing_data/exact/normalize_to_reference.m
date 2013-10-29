function [reactivity_out, reactivity_error_out] = normalize_to_reference( reactivity, reactivity_error, sequence, ref_sequence, USE_LAST, modifier );
%[reactivity_out, reactivity_error_out] = normalize_to_reference( reactivity, reactivity_error, sequence, ref_sequence, USE_LAST, modifier );
%
% Inputs:
%  reactivity       = vector of reactivity values
%  reactivity_error = vector of errors on reactivity values
%  sequence         = sequence at each reactivity position
%  ref_sequence     = [default: 'GAGUA'] sequence of referencing segment 
%  USE_LAST         = [default: 1] use last instance, rather than all instances, of reference segment.
%  modifier         = [default: ''] modifier ('DMS', 'CMCT' will trigger use of A/C, or U, respectively).
%
% (C) R. Das, Stanford University, 2013
%
if ( nargin < 4 ); help( mfilename ); return; end;
if ~exist( 'ref_sequence', 'var' ) ref_sequence = 'GAGUA'; end;
if ~exist( 'USE_LAST', 'var' ) USE_LAST = 1; end;
if ~exist( 'modifier','var') modifier = ''; end;
if length( sequence ) ~= length( reactivity ); error( 'length of sequence must match length of reactivity' ); end;

sequence = upper( sequence );
ref_sequence = upper( ref_sequence );

ref_pos_idx = strfind( sequence, ref_sequence );
if isempty( ref_pos_idx );  error( 'Could not find ref_sequence' ); end
if (USE_LAST) ref_pos_idx = ref_pos_idx(end); end;
  
ref_pos = [];
for m = 1:length( ref_pos_idx ); ref_pos = [ref_pos, ref_pos_idx(m) + [0:length(ref_sequence)-1] ]; end;

ref_pos = filter_ref_pos( ref_pos, sequence, modifier );

% only use positive...
gp = find( reactivity( ref_pos ) > 0 );
if length(gp) > 0 ;  ref_pos = ref_pos( gp ); end

[reactivity_out, norm_scalefactors, reactivity_error_out ] = ...
    quick_norm( reactivity, ref_pos, reactivity_error );  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
function ref_pos = filter_ref_pos( ref_pos, sequence, modifier );
% following should look into a database of known information about each kind of modifier. 
% future work...

if length( modifier ) == 0; return; end;
sequence = upper( sequence );
switch modifier
 case {'DMS'}
  ref_pos = filter_at_nts( ref_pos, sequence, {'A','C'} );
 case {'CMCT','CMCT (old)'}
  ref_pos = filter_at_nts( ref_pos, sequence, {'U'} );  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
function ref_pos_out = filter_at_nts( ref_pos, sequence, nts );
ref_pos_out = [];
for i = ref_pos
  if ~isempty( find( strcmp( nts, sequence( i ) ) ) ) ref_pos_out = [ref_pos_out, i ]; end;
end

