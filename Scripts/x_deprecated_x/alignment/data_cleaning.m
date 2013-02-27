function [data_all] = data_cleaning(data_all)

window = 10;
for idx = 1:size(data_all,2)
    temp_base = zeros(1, size(data_all{1},1));
    for jdx = 1000:8000
        temp_base(jdx) = min(data_all{idx}(jdx-window:jdx+window,4));
    end
    temp_base = min(-10, temp_base);
    v = find(temp_base ~= -10);
    
    data_all{idx}(v,4) = 0;
    clear temp_base;
end