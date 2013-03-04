function   [ TP_lengths, FP_lengths, TP_SHAPE, FP_SHAPE, TP_bpp, FP_bpp, TP_bpp_per_bp, FP_bpp_per_bp  ] = get_TP_FP_data( native_structure, structure, SHAPEdata, bpp, ignore_pos );
%   [ TP_lengths, FP_lengths, TP_SHAPE, FP_SHAPE, TP_bpp, FP_bpp, TP_bpp_per_bp, FP_bpp_per_bp  ] = get_TP_FP_data( native_structure, structure, SHAPEdata, bpp, ignore_pos );
%
if nargin == 0;  help( mfilename ); return; end;

TP_lengths = [];
FP_lengths = [];
TP_SHAPE   = [];
FP_SHAPE   = [];
TP_bpp = [];
FP_bpp = [];
TP_bpp_per_bp = [];
FP_bpp_per_bp = [];

native_stems = parse_stems( native_structure );

stems = parse_stems( structure );
if exist( 'ignore_pos' );  stems = filter_stems( stems, ignore_pos ); end
[ stem_correspondence, num_extra_bps_added_to_stem, num_missing_bps_in_stem ] = find_shared_stems( stems, native_stems );

FP = find( stem_correspondence == 0 );
TP = find( stem_correspondence > 0 );

stem_lengths = [];
for i = 1: length( stems );  stem_lengths(i) = length( stems{i} ); end;
TP_lengths = stem_lengths( TP );
FP_lengths = stem_lengths( FP );

SHAPE_scores = get_SHAPE_scores( stems, SHAPEdata );
TP_SHAPE = SHAPE_scores( TP );
FP_SHAPE = SHAPE_scores( FP );

stem_bpp = get_bootstrap_prob( stems, bpp );
TP_bpp = stem_bpp( TP );
FP_bpp = stem_bpp( FP );

[ TP_bpp_per_bp, FP_bpp_per_bp ] = get_TP_FP_per_bp( native_structure, structure, bpp );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function   bootstrap_prob = get_bootstrap_prob( stems, bpp );

bootstrap_prob               = zeros(  1,length( stems ) );

for j = 1:length( stems )
  stem = stems{j};
  bootprob = 0.0;

  for q = 1:size(stem,1)
    if ( bpp( stem(q,1), stem(q,2) ) > bootprob )
      bootprob = bpp( stem(q,1), stem(q,2) );
    end
  end

  bootstrap_prob(j) = bootprob;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SHAPE_scores = get_SHAPE_scores( stems, SHAPEdata );

SHAPE_scores = [];
d = SHAPEdata;

for j = 1:length( stems )
  stem = stems{j};
 
  SHAPE_score = 0.0;

  SHAPE_score = SHAPE_score + SHAPE_pseudoenergy( d(stem(1,1)) );
  SHAPE_score = SHAPE_score + SHAPE_pseudoenergy( d(stem(1,2)) );
  
  for i = 2: (size( stems,1 )-1)
    SHAPE_score = SHAPE_score + 2 * SHAPE_pseudoenergy( d(stem(i,1)) );
    SHAPE_score = SHAPE_score + 2 * SHAPE_pseudoenergy( d(stem(i,2)) );
  end
  
  SHAPE_score = SHAPE_score + SHAPE_pseudoenergy( d(stem(end,1)) );
  SHAPE_score = SHAPE_score + SHAPE_pseudoenergy( d(stem(end,2)) );

  SHAPE_scores(j) = SHAPE_score;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ TP_bpp_per_bp, FP_bpp_per_bp ] = get_TP_FP_per_bp( native_structure, structure, bpp );

TP_bpp_per_bp = [];
FP_bpp_per_bp = [];

native_structure_bps = convert_structure_to_bps( native_structure );
structure_bps =  convert_structure_to_bps( structure );
  
% could be optimized...
for q = 1:size( structure_bps, 1)
  found_match = 0;
      
  for p = 1:size( native_structure_bps, 1)
    if ( native_structure_bps( p, 1) == structure_bps( q, 1 ) & ...
	 native_structure_bps( p, 2) == structure_bps( q, 2 ) )
      found_match = 1; break;
    end
  end

  bpp_val = bpp( structure_bps( q, 1 ), structure_bps( q, 2 ) );
  if found_match
    TP_bpp_per_bp = [TP_bpp_per_bp, bpp_val ];
  else
    FP_bpp_per_bp = [FP_bpp_per_bp, bpp_val ];
  end  
end
