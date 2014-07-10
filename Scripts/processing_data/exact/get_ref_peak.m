function ref_peak = get_ref_peak( sequence, ref_segment, offset )
% ref_peak = get_ref_peak( sequence, ref_segment, offset );

if (nargin == 0 ) help( mfilename); end;
if ~exist('offset','var') || isempty(offset); offset = 0; end;

ref_pos = strfind( sequence, ref_segment );
ref_peak = [];
for i = 1:length( ref_pos );
  ref_peak = [ref_peak, ref_pos(i) - 1 + [1:length(ref_segment)] + offset ];
end