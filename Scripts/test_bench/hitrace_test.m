function hitrace_test( sample_dirs )
addpath(genpath(strcat(pwd,'/../')));

if nargin < 1
    sample_dirs = {'sample1', ...
        'sample2', ...
        'sample3', ...
        'sample4', ...
        'sample5'};
end

tic;
for i = 1:length(sample_dirs)
    a = dir(sample_dirs{i});
    nameFolds = [a(:).isdir];
    nameFolds = {a(nameFolds).name};
    
    nameFolds(ismember(nameFolds,{'.','..'})) = [];
    
    if isempty(nameFolds)
        [~, da] = quick_look(sample_dirs{i});
    else
        for j = 1:length(nameFolds)
            nameFolds{j} = [sample_dirs{i} ,'/', nameFolds{j}];
        end
        [~, da] = quick_look(nameFolds);
    end
    
    corr_test = zeros(size(da,2));
    for j = 1:size(da,2)
        for k = 1:size(da,2)
            corr_test(j,k) = corr2(da(:,j),da(:,k));
        end
    end
    
    corr_test = median(corr_test);
    
    figure;
    bar(corr_test);
    title(sprintf('sample %d - mean %f',i, mean(corr_test)));
end
toc;