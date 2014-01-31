function convert_024_030(from, to);
rdat = read_rdat_file(from);

if(strcmp(rdat.version, '0.3'))
    return;
end

structures= {};
numres = length(rdat.sequence);
for i = 1:(length(rdat.structure)/length(rdat.sequence))
    structures{i} = rdat.structure((numres * (i - 1) + 1):(numres * i));
end

if(length(structures) > 1)
    rdat.structures = structures;
end

num_structure = length(structures);
sequences = {};
str_idx = [];
for i = 1:length(rdat.data_annotations)
    data = rdat.data_annotations{i};
    idx = [];
    
    have_str = 0;
    for j = 1:length(data)
        if(strfind(upper(data{j}), 'SEQUENCE'))
            seq = strrep(upper(data{j}), 'SEQUENCE:', '');
            idx = search_Str(sequences, seq);
            if(isempty(idx))
                idx = length(sequences) + 1;
                sequences{idx} = seq;                
            end
            data{j} = sprintf('sequence:%d',idx);
        elseif(strfind(upper(data{j}), 'MODIFIER'))
            modifier = strrep(upper(data{j}), 'MODIFIER:', '');
            if((strcmp(modifier, 'SHAPE') || strcmp(modifier, 'DMS') || strcmp(modifier,'CMCT')))
                have_str = 1;
            end
        end
    end
    
    if(~isempty(idx))
        if(length(str_idx) < idx)
            str_idx(idx) = have_str;
        else
            str_idx(idx) = str_idx(idx) + have_str;
        end
        str = mod(str_idx(idx), num_structure);
        if(str == 0)
            str = num_structure;
        end
        data{length(data)+1} = sprintf('structure:%d',str);        
    end
    
    rdat.data_annotations{i} = data;
end

rdat.sequences = sequences;
rdat.version = '0.3';
output_rdat03_to_file(to, rdat);

function idx = search_Str(strs, str)
idx = [];
for i = 1:length(strs)
    if(strcmp(strs{i}, str))
        idx = [idx, i];
    end
end