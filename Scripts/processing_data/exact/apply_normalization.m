function [reactivity, reactivity_error, norm_scalefactors ] = apply_normalization( reactionProb, reactionProb_error, refpos, seqpos, data_type, sequence, offset );
%
%  [reactivity, reactivity_error ] = apply_normalization( reactionProb, reactionProb_error, refpos, seqpos, data_type, sequence, offset );
%
% Data type of DMS will flag use of A & C; CMCT will flag use of U.
% Data type of 'nomod' or ddATP, ddCTP, ddGTP, or ddTTP will flag use of
% a  verage normalization scalefactor used in other lanes.
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
%  norm_scalefactors   = scalefactors applied to each lane
%
% (C) R. Das, Stanford University, 2013, 2018

for j = 1:size( reactionProb, 2 )
    
    ref_pos_in = [];
    for k = 1:length( refpos )
        if ~isempty( data_type ) & ~isempty( sequence );
            nt = upper(sequence( refpos(k) - offset ));
            if ( strcmp( data_type{j}, 'DMS' )  & nt ~= 'A' & nt ~= 'C' ); continue; end;
            if ( strcmp( data_type{j}, 'CMCT' ) & nt ~= 'U' ); continue; end;
           if ( strcmp( data_type{j}, 'glyoxal' ) & nt ~= 'G' ); continue; end;
        end
        ref_pos_in = [ ref_pos_in, find( seqpos == refpos(k) ) ];
    end
    
    [reactivity(:,j), norm_scalefactors(j), reactivity_error(:,j) ] = quick_norm( reactionProb(:,j), ref_pos_in, reactionProb_error(:,j) );
    
end

% don't normalize the 'reference' lanes ? nomod, ddNTP ladders ? instead
% apply the normalization factor used in other lanes.
reference_lanes = find( strcmp( data_type, 'nomod' ) | strcmp( data_type, 'ddATP' )| strcmp( data_type, 'ddTTP' ) | strcmp( data_type, 'ddCTP' ) |  strcmp( data_type, 'ddGTP' ) );
normalize_lanes = setdiff( [1:size( reactionProb, 2 )], reference_lanes );
mean_scalefactor = 0;
if length( normalize_lanes ) > 0; mean_scalefactor = mean( norm_scalefactors( normalize_lanes ) ); end;
for j = reference_lanes
    norm_scalefactors(j) = mean_scalefactor;
    reactivity(:,j) = reactionProb(:,j) * mean_scalefactor;
    reactivity_error(:,j) = reactionProb_error(:,j) * mean_scalefactor;
end
