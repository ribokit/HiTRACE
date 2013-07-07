function hitrace_test()
addpath(genpath(strcat(pwd,'/../')));

sample_dirs = {{'sample1'}, ...
    {'sample2/1_data', 'sample2/2_nomod_ladder'}, ...
    {'sample3'}, ...
    {'sample4'}};

sameple_seq = { 'GGAAAGCUGACAGGAUAUGGGCAUGACAAAAGUCAUGCGCUGGCUACAAAAGUAGGCAGCCUAGAUCAAAAGAUCUAGCCAGAAGGGUCAGCAAAGAAACAACAACAACAAC', ...
    'GGCAGGACGAUGUUGCUAUGAAUGUAUAGUACGACAGGAUAUUUGUUAUUAAAGAAGGGUCGCCGCGUACACCUAUGCGGACAUCGUUAAAGAAACAACAACAACAAC', ...
    'GGCGCGUUUGAAGGAUAUCGGCCAGCGGCAUCGUUGGUAACAUGUAGUCGGCUAUUUGUUAACUGUUCUAUAACGGUUCGAGAAGGUCAAACAAAGAAACAACAACAACAAC',...
    'GGAAAGUAGAGGAUAUGCAGUUCAGUCAAAAGACAGAGAAAACUCAGAACAGUACACAGAAAACUGAGUGAAAACACAGUACAGCAGAAGGCUACAAAGAAACAACAACAACAAC'};

sample_str = { '.....((((((......((((((((((....)))))))((((.((((....)))).))))(((((((....)))))))))).....))))))....................', ...
    '.....(((((((((((((((....)))))))((((......((((....)))).....))))(((((((....)))))))))))))))....................', ...
    '.....((((((......((((((((((....)))))))((((.((((....)))).))))(((((((....)))))))))).....))))))....................', ...
    '.....((((......(((.((((.(((....))).(((....))).)))).((((.(((....))).(((....))).)))).))).....))))....................'};

sample_type = {'SHAPE', 'SHAPE', 'DMS', 'DMS', 'nomod', 'ddTTP'};

tic;
for i = 1:length(sample_dirs)
    [d{i}, da{i}] = quick_look(sample_dirs{i});
    
    corr_test = zeros(size(da{i},2));
    for j = 1:size(da{i},2)
        for k = 1:size(da{i},2)
            corr_test(j,k) = corr2(da{i}(:,j),da{i}(:,k));
        end
    end
    
    corr_test = median(corr_test);
    
    figure;
    bar(corr_test);
    title(sprintf('sample %d - mean %f',i, mean(corr_test)));
    
    
end
toc;
