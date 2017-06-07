function mutpos = generate_mutpos_from_rdat(cell_array,rdat)
% mutpos = GENERATE_MUTPOS_FROM_RDAT(data_annotations)
%
% Generates mutpos array from {data_annotations} from rdat file.
% WT will be recorded as NaN.
%
% by T47, Apr 2013
% replaced by same function used in read_rdat_file() in RDATkit
% (get_mutation_info_from_tag) by Rhiju, 2017.
%

if nargin == 0; help( mfilename ); return; end;
mutation_tags = get_tag( cell_array, 'mutation' );
mutpos = ones(1, size(cell_array, 2)) * NaN;
for i = 2:size(mutpos, 2)
    [mutpos_inferred,~] = get_mutation_info_from_tag( mutation_tags{i}, rdat );
    if ~isempty( mutpos_inferred )
        mutpos(i) = mutpos_inferred;
    end
end;
%fprintf(['mutpos (1 x ', num2str(length(mutpos)), ') is generated from d_rdat.data_annotations.\n']);
