function image_diagram(imagex)

%
% IMAGE_DIAGRAM(imagex);
%
% Shows imagex in a full-screened equal-axis figure.
%
% by T47, May 2013.
%

hold off; clf;
figure_full_screen;
image(imagex / 256); 
hold on; axis equal; axis off;