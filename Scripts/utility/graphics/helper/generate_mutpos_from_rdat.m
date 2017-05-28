function mutpos = generate_mutpos_from_rdat(cell_array)

% mutpos = GENERATE_MUTPOS_FROM_RDAT(data_annotations)
%
% Generates mutpos array from {data_annotations} from rdat file.
% WT will be recorded as NaN.
%
% by T47, Apr 2013
%

if nargin == 0; help( mfilename ); return; end;

mutpos = ones(1, size(cell_array, 2)) * NaN; 
for i = 2:size(mutpos, 2)
  mutpos_inferred = str2num(strrep(strrep(strrep(strrep(strrep(strrep(cell_array{i}{1}, 'mutation:', ''), 'G', ''), 'A', ''), 'C', ''), 'U', ''), 'X', ''));
  if ~isempty( mutpos_inferred )
    mutpos(i) = mutpos_inferred;
  end
end;
fprintf(['mutpos (1 x ', num2str(length(mutpos)), ') is generated from d_rdat.data_annotations.\n']);

% mutpos