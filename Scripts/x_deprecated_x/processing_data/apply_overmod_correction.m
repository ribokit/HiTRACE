function  [ d_unmod, correction ] = apply_overmod_correction( d, modification_ratio );
% APPLY_OVERMOD_CORRECTION: applies approximate overmodification correction
%
%  [ d_unmod, correction ] = apply_overmod_correction( d, modification_ratio );
%
% Deprecated soon.
%
if nargin == 0;  help( mfilename ); return; end;


% Don't read the first few residues -- dominated by fluorescence
% from "excess" primers that were unextended.
eat_in = 4;
d_cumulative_sum = 0 * d;

% d_cumulative_sum integrates band intensity from smaller to larger DNAs.
d_cumulative_sum( eat_in:end ) = cumsum( d(eat_in:end) );
d_cumulative_sum = d_cumulative_sum/d_cumulative_sum( end );

%plot( d_cumulative_sum ); pause;
correction = exp ( modification_ratio .* d_cumulative_sum );
d_unmod = d .* correction;
