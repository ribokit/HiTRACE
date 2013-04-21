function print_bpp_Z (bpp, Z, scale_factor, filename, if_print)

% PRINT_BPP_Z(bpp, Z, [scale_factor], [filename], [if_print])
%
% Prints the base pairing probability matrix and Z score matrix in 600*600 squared figures, and
%    output .png images if asked.
%
% Arguments
% =========
%   bpp                  Required        Provides the base pairing probability matrix.
%   Z                    Required        Provides the Z score matrix.
%   [scale_factor]       Optional        Provides the scale factor for Z score matrix,  
%                                            value should be negative. Default is -5.
%   [filename]           Optional        Provides the filename for output images. Default
%                                            is blank, i.e. 'Z.png' and 'bpp.png'. Prefix, 
%                                            underlines, and '.png' extension automatically
%                                            append.
%   [if_print]           Optional        Provides whether to output image file or not.
%                                            Default is 1 (YES).
%
% Notes
% =====
% All current opened figures will be lost. Save before run.
%
% by T47, Apr 2013
%

if ~exist('scale_factor') || isempty('scale_factor') || scale_factor == 0; scale_factor = -5; end;
if ~exist('filename') || isempty('filename'); filename = ''; else; filename = ['_', filename]; end;
if ~exist('if_print') || isempty('if_print'); if_print = 1; end;
fprintf('scale_factor = %d\n', scale_factor);

close all;

% plot Z score
h1 = figure(1);
set(h1,'Position',[0, 0, 600, 600]);
set(gcf, 'PaperOrientation', 'landscape', 'PaperPositionMode', 'auto', 'color', 'white');
image(Z * scale_factor); colormap(1-gray());
if if_print == 1; print(h1,'-dpng',['Z', filename, '.png']); end;

% plot bpp
h2 = figure(2);
image(bpp*100); colormap(1-gray());
set(h2,'Position',[100, 0, 600, 600]);
set(gcf, 'PaperOrientation', 'landscape', 'PaperPositionMode', 'auto', 'color', 'white');
image(bpp*100); colormap(1-gray());
if if_print == 1; print(h2,'-dpng',['bpp', filename, '.png']); end;



