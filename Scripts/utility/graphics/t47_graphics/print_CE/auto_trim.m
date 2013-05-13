function [d_trim, xsel_trim] = auto_trim(d_org, xsel_org, up_offset, low_offset, flg)

d_upper_bound = max([round((min(xsel_org) - up_offset)), 1]);
d_lower_bound = min([round((max(xsel_org) + low_offset)), size(d_org, 1)]);

if flg == 1;
    d_trim = d_org(d_upper_bound:d_lower_bound, :);
    xsel_trim = xsel_org - d_upper_bound + 1;
    
    fprintf([', with upper offset (', num2str(up_offset), ') and lower offset (', ...
        num2str(low_offset), '), and trimmed to ', num2str(d_upper_bound), ' : ', num2str(d_lower_bound)]);
end;
