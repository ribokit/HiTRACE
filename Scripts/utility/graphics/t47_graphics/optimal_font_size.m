function ft_sz = optimal_font_size(text_obj, max_l, max_h)
set(text_obj, 'Units', 'Inches')
[text_l, text_h, text_fs] = get_text_size(text_obj);

while (text_l < max_l) && (text_h < max_h)
    text_fs = text_fs + 0.5;
    set(text_obj, 'FontSize', text_fs);
    [text_l, text_h, text_fs] = get_text_size(text_obj);
end;

while (text_l > max_l) && (text_h > max_h)
    text_fs = text_fs - 0.5;
    set(text_obj, 'FontSize', text_fs);
    [text_l, text_h, text_fs] = get_text_size(text_obj);
end;

set(text_obj, 'FontSize', text_fs);
ft_sz = text_fs;

%%%%
function [t_l, t_h, t_fs] = get_text_size(t_obj)

t_xt = get(t_obj, 'Extent');
t_l = t_xt(3);
t_h = t_xt(4);
t_fs = get(t_obj, 'FontSize');
