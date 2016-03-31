function show_boxes(imagex, residue_locations, square_width, is_round)
% show_boxes(imagex, residue_locations, square_width, is_round)
%
% imagex            = RGB image [M1 x M2 x 3 matrix] read in from, say a tif file with
%                      the 'imread' command.
% residue_locations = 2 x N matrix with the (x,y) positions of each 'sequence' position
%                      on the image. Initially you can set this to the empty set [].
% square_width      = [default 24] how big to make squares.
% is_round          = [default 0] whether to use squares or circles.
%
%

figure(1); subplot(1,1,1); hold off; image(imagex); hold on;
[xsize, ysize, zsize] = size(imagex);
axis([0 ysize 0 xsize]); zoomedin = 0;
figure_full_screen(); axis equal;

if ~exist('square_width', 'var') || isempty(square_width); square_width = 40; end;
if ~exist('is_round', 'var') || isempty(is_round); is_round = 0; end;

count = length(residue_locations) + 1;
for k = 1:(count - 1)
    xpick = residue_locations(1, k);
    ypick = residue_locations(2, k);
    if is_round;
        h(k) = rectangle('Position', [xpick - square_width/2, ypick-square_width/2, square_width, square_width], ...
            'Curvature', [1 1]);
    else
        h(k) = rectangle('Position', [xpick - square_width/2, ypick-square_width/2, square_width, square_width]);
    end;
    set(h(k), 'edgecolor', 'b');
end
