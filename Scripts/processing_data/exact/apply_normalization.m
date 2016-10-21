function [reactivity, reactivity_error ] = apply_normalization( reactionProb, reactionProb_error, refpos, seqpos, data_type, sequence, offset );
%
%  [reactivity, reactivity_error ] = apply_normalization( reactionProb, reactionProb_error, refpos, seqpos, data_type, sequence, offset );
%
% Required Inputs
%  reactionProb        = reactivities
%  reactionProb_error  = errors on reactivities
%  refpos              = reference positions; after normalization these will average to 1.0.
%
% Optional Inputs:
%  seqpos              = conventional sequence numbers for input data (default: 1, 2, ...)
%  data_type           = cell of strings {'DMS','CMCT',...} -- DMS will flag use of A & C; CMCT will flag use of U
%  sequence            = RNA sequence
%  offset              = offset to add to 1, 2, ... to conventional numbering in sequence.
%
% Outputs
%  reactivity          = normalized reactivity
%  reactivity_error    = error on normalized reactivity
%
% (C) R. Das, Stanford University, 2013

for j = 1:size( reactionProb, 2 )
    
    ref_pos_in = [];
    for k = 1:length( refpos )
        if ~isempty( data_type ) & ~isempty( sequence );
            nt = upper(sequence( refpos(k) - offset ));
            if ( strcmp( data_type{j}, 'DMS' )  & nt ~= 'A' & nt ~= 'C' ); continue; end;
            if ( strcmp( data_type{j}, 'CMCT' ) & nt ~= 'U' ); continue; end;
        end
        ref_pos_in = [ ref_pos_in, find( seqpos == refpos(k) ) ];
    end
    
    [reactivity(:,j), norm_scalefactors, reactivity_error(:,j) ] = quick_norm( reactionProb(:,j), ref_pos_in, reactionProb_error(:,j) );
    
end