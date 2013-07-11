function hitrace_test()
addpath(genpath(strcat(pwd,'/../')));

load('./compare/original.mat');

sample_dirs = {{'sample1'}, ...
    {'sample2/1_data', 'sample2/2_nomod_ladder'}, ...
    {'sample3'}, ...
    {'sample4'}};

sample_seq = { 'GGAAAGCUGACAGGAUAUGGGCAUGACAAAAGUCAUGCGCUGGCUACAAAAGUAGGCAGCCUAGAUCAAAAGAUCUAGCCAGAAGGGUCAGCAAAGAAACAACAACAACAAC', ...
    'GGCAGGACGAUGUUGCUAUGAAUGUAUAGUACGACAGGAUAUUUGUUAUUAAAGAAGGGUCGCCGCGUACACCUAUGCGGACAUCGUUAAAGAAACAACAACAACAAC', ...
    'GGCGCGUUUGAAGGAUAUCGGCCAGCGGCAUCGUUGGUAACAUGUAGUCGGCUAUUUGUUAACUGUUCUAUAACGGUUCGAGAAGGUCAAACAAAGAAACAACAACAACAAC',...
    'GGAAAGUAGAGGAUAUGCAGUUCAGUCAAAAGACAGAGAAAACUCAGAACAGUACACAGAAAACUGAGUGAAAACACAGUACAGCAGAAGGCUACAAAGAAACAACAACAACAAC'};

sample_str = { '.....((((((......((((((((((....)))))))((((.((((....)))).))))(((((((....)))))))))).....))))))....................', ...
    '.....(((((((((((((((....)))))))((((......((((....)))).....))))(((((((....)))))))))))))))....................', ...
    '.....((((((......((((((((((....)))))))((((.((((....)))).))))(((((((....)))))))))).....))))))....................', ...
    '.....((((......(((.((((.(((....))).(((....))).)))).((((.(((....))).(((....))).)))).))).....))))....................'};

sample_type = {{'SHAPE', 'SHAPE', 'DMS', 'DMS', 'nomod', 'ddTTP'},...
    {'DMS', 'DMS', 'DMS', 'nomod', 'ddTTP'},...
    {'SHAPE', 'SHAPE', 'DMS', 'DMS', 'nomod', 'ddTTP'},...
    {'SHAPE', 'SHAPE', 'DMS', 'DMS', 'nomod', 'ddTTP'}};

tic;
for i = 1:length(sample_dirs)
    [d{i}, da{i}] = quick_look(sample_dirs{i});
    
    corr_test = zeros(size(da{i},2));
    corr_test_org = zeros(size(da_org{i},2));
    for j = 1:size(da{i},2)
        for k = 1:size(da{i},2)
            corr_test(j,k) = corr2(da{i}(:,j),da{i}(:,k));
            corr_test_org(j,k) = corr2(da_org{i}(:,j),da_org{i}(:,k));
        end
    end
    
    corr_test = median(corr_test);
    corr_test_org = median(corr_test_org);
    
    figure;
    bar(corr_test);
    title(sprintf('sample %d - mean %f',i, mean(corr_test)));   
    
    corr_mean(i) = mean(corr_test);
    corr_mean_org(i) = mean(corr_test_org);
    
    area_pred{i} = generate_area_pred(sample_seq{i}(1:end-20), sample_str{i}(1:end-20), 0, sample_type{i}, length(sample_type{i}));
    
    area_pred_reverse = area_pred{i}(end:-1:1,:);
    
    xsel{i} = auto_assign_sequence(d{i}, sample_seq{i}(1:end-20), area_pred_reverse, 0, [], sample_type{i});    
    area_peak{i} = fit_to_gaussians( d{i}, xsel{i} );
    
    corr_area = zeros(1,size(da{i},2));
    for j = 1:size(da{i},2)
        corr_area(j) = corr2(area_peak{i}(:,j),area_peak_org{i}(:,j));
    end
    
    corr_area_mean(i) = mean(corr_area);
end
toc;

fprintf(1,'Report... \n\n');
for i = 1:length(sample_dirs)
    fprintf(1,'Sample %d Align Correlation Coeff. - current code: %f, original code: %f\n', i, corr_mean(i), corr_mean_org(i)); 
    fprintf(1,'Sample %d Fitting Area Correlation Coeff. to original code - %f\n\n', i, corr_area_mean(i)); 
end