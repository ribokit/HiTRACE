function annotation_subset = annotation_type_finder(annotations, type_string)

annotation_subset = '';
for i = 1:length(annotations)
    tok = strfind(annotations{i}, type_string);
    if ~isempty(tok);
        annotation_subset = [annotation_subset, '  ', ...
            annotations{i}((tok + length(type_string) + 1):end)];
    end;
end;

if ~isempty(annotation_subset); 
    annotation_subset = annotation_subset(3:end);
end;