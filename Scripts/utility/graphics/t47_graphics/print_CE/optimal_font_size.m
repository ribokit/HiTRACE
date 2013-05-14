function ft_sz = optimal_font_size(text_obj, max_l, max_h)

% ft_sz = OPTIMAL_FONT_SIZE(text_obj, max_l, max_h);
%
% Find optimal text font size for a certain text object under defined
%   maximum height and width.
%
% Input
% =====
%   text_obj    Required        Provides the text object handle. Feed in
%                               	the handle instead of string.
%   max_l       Required        Provides the maximum width of text field.
%   max_h       Required        Provides the maximum height of text field.
%
% Output
% ======
%   ft_sz                       Provides the optimal fit font size.
%
% Note
% ====
% Feed 0 for max_l or max_h to disable the constraint.
% Text label will be adjust to ft_sz automatically on figure.
%
% by T47, May 2013.
%

if nargin == 0; help( mfilename ); return; end;

if max_l == 0 && max_h == 0;
    ft_sz = 0;
    return;
end;

set(text_obj, 'Units', 'Inches')
[~, ~, text_fs] = get_text_size(text_obj);

% increase ft_sz if smaller than maximum
condition = get_condition(text_obj, max_l, max_h, -1);
while condition
    text_fs = text_fs + 0.5;
    set(text_obj, 'FontSize', text_fs);
    condition = get_condition(text_obj, max_l, max_h, -1);
end;

% decrease ft_sz if bigger than maximum
condition = get_condition(text_obj, max_l, max_h, 1);
while condition
    text_fs = text_fs - 0.5;
    set(text_obj, 'FontSize', text_fs);
    condition = get_condition(text_obj, max_l, max_h, 1);
end;

set(text_obj, 'FontSize', text_fs);
ft_sz = text_fs;

%%%%
% get text width, height and font size from handle
function [t_l, t_h, t_fs] = get_text_size(t_obj)

t_xt = get(t_obj, 'Extent');
t_l = t_xt(3);
t_h = t_xt(4);
t_fs = get(t_obj, 'FontSize');

%%%%
% condtion of while loop; -1 for if smaller, +1 for if bigger
function cdt = get_condition(t_obj, m_l, m_h, drct)

[t_l, t_h] = get_text_size(t_obj);

if drct == -1;
    if m_l ~= 0 && m_h ~= 0 ;
        cdt = (t_l < m_l) && (t_h < m_h);
    elseif m_l ~= 0
        cdt = t_l < m_l;
    elseif m_h ~= 0
        cdt = t_h < m_h;
    end;
elseif drct == 1;
    if m_l ~= 0 && m_h ~= 0 ;
        cdt = (t_l > m_l) && (t_h > m_h);
    elseif m_l ~= 0
        cdt = t_l > m_l;
    elseif m_h ~= 0
        cdt = t_h > m_h;
    end;
end;    
