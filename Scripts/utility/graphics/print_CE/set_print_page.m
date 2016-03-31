function set_print_page(fig, orient, pos, txt)

%
% SET_PRINT_PAGE(fig, orient, pos, txt);
%
% Set figure to optimal printing size on letter papers.
%
% Input
% =====
%   fig         Required    Provides the handle for figure, e.g. figure(1)
%                               or gcf.
%   orient      Optional    Provides the orientation of pages. 0 for
%                               landscape, 1 for portrait. Default is 1.
%   pos         Optional    Provides 1x4 double array for figure window
%                               position and size on screen. Default is 
%                               [0, 0, 600, 800];
%   txt         Optional    Provides text string for figure window title.
%
% 
% by T47, May 2013.
%

if nargin == 0; help( mfilename ); return; end;

if exist('txt','var'); set(fig, 'Name', txt); end;
if ~exist('pos','var') || isempty(pos); pos = [0 0 600 800]; end;

set(fig, 'Position', pos);
set(fig, 'PaperPositionMode', 'Manual', 'Color', 'White');

if ~exist('orient','var') || isempty(orient); orient = 1; end;
orient = is_valid_boolean(orient);
if orient == 0;
    set(fig, 'PaperOrientation', 'Landscape', ...
        'PaperSize', [11 8.5], 'PaperPosition', [-0.65 0.15 12 8]);
else
    set(fig, 'PaperOrientation', 'Portrait', ...
        'PaperSize', [8.5 11], 'PaperPosition', [0 1 8.5 10.5]);
end;

