function [pos_sub, txt_sub] = parse_xsel_label(xsel_cell, h_sub, ct_offset)

pos_sub = []; txt_sub = {}; 

for k = 1:size(xsel_cell, 2)
    if ~isempty(xsel_cell{h_sub, k, 2});
        pos_sub(ct_offset) = xsel_cell{h_sub, k, 2};
        txt_sub{ct_offset} = xsel_cell{h_sub, k, 1};
        ct_offset = ct_offset + 1;
    end;
end;
