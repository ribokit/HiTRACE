function im_txt = rasterize_text(txt, ft_sz, clr, ft)

%
% text_image = RASTERIZE_TEXT(text, font_size, font_color, font_name);
%
% Returns rasterized text image (N x N x 3 unit8] of defined text with size
%  and color.
% A window will pop out and disappear during running, which is annoying.
%
%
% Input
% =====
%   text            Required        Provides the text string.
%   font_size       Required        Provides the font size number.
%   font_color      Optional        Provides the font color code. Default
%                                    is black [0 0 0].
%   font_name       Optional        Provides the font name. Default is
%                                    'Helvetica'.
%
% Output
% ======
%   text_image                      Gives the image matrix of text.
%
%
% by T47, May 2013.
%

if nargin == 0; help( mfilename ); return; end;

if ~exist('clr','var') || isempty(clr); clr = 'k'; end;
if ~exist('ft','var') || isempty(ft); ft = 'Helvetica'; end;

h = figure;
set(gca, 'XTickLabel', '', 'YTickLabel', '');

% scale up font size to overcome axis size problem
l_offset = 1 + 1 / length(txt);
if ft_sz <= 60 && length(txt) <= 4; scale_factor = 10; else; scale_factor = 5; end;

h_txt = text(0, 0, txt, ...
    'HorizontalAlignment', 'Left', 'VerticalAlignment', 'Bottom', ...
    'FontSize', ft_sz * scale_factor, 'Color', clr, ...
    'BackgroundColor', 'White');
set(h_txt, 'Units', 'Pixels', 'FontName', ft);

% set figure size to minimum and take a snapshot
h_ext = get(h_txt, 'Extent');
h_l = h_ext(3); h_h = h_ext(4);
set(h, 'Color', 'White', 'Units', 'Pixels', ...
    'Position', [0, 0, h_l * l_offset, h_h]);

frame_txt = getframe(gcf);
close(h);
image_txt = frame2im(frame_txt);

% resize back to 1x
im_txt = imresize(image_txt, 1 / scale_factor);
