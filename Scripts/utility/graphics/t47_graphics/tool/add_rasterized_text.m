function imagex_text = add_rasterized_text(imagex, pos, txt, ft_sz, clr, ft)

%
% imagex_test = ADD_RASTERIZED_TEXT(imagex, position, text, ...
%                                   font_size, font_color, font_name);
%
% Returns rasterized text image (N x N x 3 unit8] of defined text with size
%  and color added to input image at defined position.
% A window will pop out and disappear during running, which is annoying.
%
%
% Input
% =====
%   imagex          Required        Provides the original image matrix.
%   position        Required        Provides the coordinate to add the text
%                                    image into original image.
%   text            Required        Provides the text string.
%   font_size       Required        Provides the font size number.
%   font_color      Optional        Provides the font color code. Default
%                                    is black [0 0 0].
%   font_name       Optional        Provides the font name. Default is
%                                    'Helvetica'.
%
% Output
% ======
%   imagex_text                     Gives the image matrix with text added.
%
%
% by T47, May 2013.
%

if nargin == 0; help( mfilename ); return; end;

if ~exist('clr','var') || isempty(clr); clr = 'k'; end;
if ~exist('ft','var') || isempty(ft); ft = 'Helvetica'; end;

imagex_text = imagex;

% get text image
lbl = double(rasterize_text(txt, ft_sz, clr, ft));

% add to position
xbins = (pos(2) + 1):(pos(2) + size(lbl, 1));
ybins = (pos(1) + 1):(pos(1) + size(lbl, 2));
imagex_text(xbins, ybins, :) = lbl;
