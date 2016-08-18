function figure_full_screen(fig)

%
% FIGURE_FULL_SCREEN(fig);
%
% Sets figure to full screen display. Feed in figure handle, e.g. figure(1)
% or gcf. Default is gcf.
%
% by T47, May 2013.
%

if ~exist('fig','var') || isempty(fig); fig = gcf; end;

set(fig, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
