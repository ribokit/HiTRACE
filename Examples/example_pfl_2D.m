%% STEP 1
filenames = {'pfl_2D_SHAPE_plus_409193'};
[d_align, d_ref_align, ylimit, labels] = quick_look(filenames, [], 1:72);

%% STEP 2
% excludes the bad lane #56
align_blocks = {[1:55, 57:72]};
d_align_before_more_alignment = d_align;
d_align_dp_fine = align_by_DP_fine(d_align_before_more_alignment, align_blocks);
d_align = d_align_dp_fine;


%% STEP 3
sequence = 'GGAGACCTCGAGTAGAGGTCAAAAGGGTCGTGACTGGCGAACAGGTGGGAAACCACCGGGGAGCGACCCCGGCATCGATAGCCGCCCGCCTGGGCAAACAACTCGAGTAGAGTTGACAACAAAGAAACAACAACAACAAC';
sequence = strrep(sequence, 'T', 'U');
% use crystallographic struture for guidance
% untested ones should not be used here
structure = '...((((((.....))))))....(((((((....[[[[....(((((....))))).....)))))))...........(((..]]]]...)))...((((((.....)))))).........................';
offset = -24;
first_RT_nucleotide = length(sequence) - 20 + offset; 

% the first lane is WT, thus NaN in terms of mutation
data_types{1} = 'NaN';
for i = 2:72;
    data_types{i} = [num2str(i - 1)];
end;


%% STEP 4
xsel = []; clf;
% remove the extra final band
[xsel, seqpos, area_pred] = annotate_sequence(d_align, xsel, sequence, offset, data_types, first_RT_nucleotide, structure);

%% STEP 5
[area_peak, darea_peak] = fit_to_gaussians(d_align, xsel);
area_peak(:, 56) = 0;
darea_peak(:, 56) = 0;


%% STEP 9
filename = 'pfl_2D_SHAPE_plus.rdat';
name = 'RNA Puzzle 13';
comments = { ...
    'Preliminary data', ...
    'RNA was heat to 90 C for 2 min and cool on ice for 2 min, then incubated at 37 C for 20 min at pH 8.0,', ...
    '  in presence of 10 uM ZMP ligand with 10 mM MgCl2 to aid folding.', ...
    'Note: ', ...
    'Taken in early May 2015 for thirteenth RNA Puzzle community-wide modeling trial.', ...
};

annotations = { ...
    'experimentType:MutateAndMap', ...
    'chemical:Na-HEPES:50mM(pH8.0)', 'chemical:MgCl2:10mM', 'chemical:ZMP:10uM', ...
    'temperature:24C', ...
    'modifier:SHAPE', ...
};

% loop generates a cell of mutation labels, e.g. U55A
data_annotations{1} = 'mutation:WT';
for j = 2:size(area_peak, 2); 
    data_annotations{j} = {['mutation:', sequence(j - 1 - offset), num2str(j - 1), DNA2RNA(complement(sequence(j - 1 - offset)))]};
end;
% manually mark bad lane
data_annotations{56} = {'mutation:U55A', 'warning:badQuality'};

output_workspace_to_rdat_file(filename, name, sequence, offset, seqpos, area_peak, ...
    structure, annotations, data_annotations, darea_peak, [], [], [], comments);
d_rdat = show_rdat(filename);


%% STEP 8
d_align(:, 56) = 0;
% you need the RDAT instance for printing
print_CE_split(d_align, d_rdat, [], [], xsel, [], [], area_pred, ...
    [2 2], [1 0 1 1 1 1 1 0 1 1 1], ...
    {'', 'print_CE_split_output', 'SHAPE', 'T47', 'May 2015'});


%% STEP 10
Z = output_Zscore_from_rdat('pfl_SHAPE_2Dbonus.txt', {filename});

% filtering of Z score by sequence positions
Z_cutoff_mean = Z;
for i = 1:size(Z_cutoff_mean, 1);
    Z_cutoff_mean(i, Z_cutoff_mean(i, :) >= mean(Z_cutoff_mean(i, :))) = 0;
end;

Z_cutoff_1std = Z;
for i = 1:size(Z_cutoff_1std, 1);
    Z_cutoff_1std(i, Z_cutoff_1std(i, :) >= -std(Z_cutoff_1std(i, :)) + mean(Z_cutoff_1std(i, :))) = 0;
end;


%% STEP 11
% run Fold with 100x bootstrap on 2D SHAPE data
[structure_2D_Fold_SHAPE, bpp_2D_Fold_SHAPE] = rna_structure(sequence, [], offset, seqpos, Z, 100, 0);
% run ShapeKnot with 100x bootstrap on 2D SHAPE data
[structure_2D_Spkt_SHAPE, bpp_2D_Spkt_SHAPE] = rna_structure(sequence, [], offset, seqpos, Z, 100, 1);

%% STEP 12
output_varna_html('pfl_2D_Fold_SHAPE.html', sequence, structure_2D_Fold_SHAPE, structure, structure_2D_Fold_SHAPE, offset, [], [], [], bpp_2D_Fold_SHAPE);
output_varna_html('pfl_2D_Spkt_SHAPE.html', sequence, structure_2D_Spkt_SHAPE, structure, structure_2D_Spkt_SHAPE, offset, [], [], [], bpp_2D_Spkt_SHAPE);

print_bpp_Z(bpp_2D_Fold_SHAPE, Z, -15, '2D_Fold_SHAPE');
print_bpp_Z(bpp_2D_Spkt_SHAPE, Z, -15, '2D_Spkt_SHAPE');

