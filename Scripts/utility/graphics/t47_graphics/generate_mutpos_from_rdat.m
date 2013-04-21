function mutpos = generate_mutpos_from_rdat(cell_array)

% mutpos = GENERATE_MUTPOS_FROM_RDAT(data_annotations)
%
% Generates mutpos array from {data_annotations} from rdat file.
% WT will be recorded as NaN.
%
% by T47, Apr 2013
%

mutpos = zeros(1, size(cell_array, 2)); 
mutpos(1) = NaN;
for i = 2:size(mutpos, 2)
    mutpos(i) = str2num(strrep(strrep(strrep(strrep(strrep(cell_array{i}{1}, 'mutation:', ''), 'G', ''), 'A', ''), 'C', ''), 'U', ''));
end;
fprintf('\n');
fprintf(['mutpos (1 x ', num2str(length(mutpos)), ') is generated from d_rdat.data_annotations.\n']);