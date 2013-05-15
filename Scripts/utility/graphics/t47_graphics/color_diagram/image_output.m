function image_output(imagex, file_name, res)

%
% IMAGE_OUTPUT(imagex, file_name, resolution)
%
% Output colored diagram to hi-resolution tiff image file.
% 
%
% Input
% =====
% imagex            = RGB image [M1 x M2 x 3 matrix] read in from, say a tif 
%                      file with the 'imread' command.
% file_name         = [default 'color_diagram_output.tiff'] file name for
%                      image file. Extension will automatically append.
% resolution        = [default 300] resolution of tiff.
%
%
% by T47, May 2013.
%

if nargin == 0; help( mfilename ); return; end;

if ~exist('file_name','var') || isempty(file_name); file_name = 'color_diagram_output'; end;
file_name = [file_name, '.tiff'];
if ~exist('res','var') || isempty(res); res = 300; end;

% output to hi-res non-compressed tiff
d_tiff = im2uint8(imagex / 256);
imwrite(d_tiff, file_name, 'TIFF', 'Resolution', res, 'Compression', 'lzw');

