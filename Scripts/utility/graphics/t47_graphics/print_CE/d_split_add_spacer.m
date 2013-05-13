function d_spl = d_split_add_spacer(d_org, H, h_length, h_sp, W, w_length, w_sp)

d_spl = zeros([W, H, (h_length + 2 * h_sp), (w_length + 2 * w_sp)]); 
for i = 1:W
    for j = 1:H
        d_spl(i, j, ((h_sp + 1):(h_sp + h_length)), ((w_sp + 1):(w_sp + w_length))) = ...
            d_org((h_length * (j - 1) + 1):(h_length * j),...
            (w_length * (i - 1) + 1):(w_length * i));
    end;
end;
