function make_lines(line_pos,color_code,line_width, low_bd, high_bd)
%
% MAKE_LINES([line_pos], [color_code], [line_width], [low_bd], [high_bd]);
%
% Draws vertical lines in a figure.
%
% Input
% =====
%   line_pos        Denotes the position of vertical lines. Format in
%                       double array. Default is [-?:+?] on x-axis.
%   color_code      Denotes color of lines. Default is 'k' (black).
%   line_width      Denotes width of lines. Default is 1.
%   low_bd          Denotes whether draw lines to -? on y-axis. Default is
%                       1 (YES).
%   high_bd         Denotes whether draw lines to +? on y-axis. Default is
%                       1 (YES).
%
% e.g. MAKE_LINES() will draw black vertical lines on every integer number
%       on the figure from bottom to top with line width 1.
%
% by Rhiju Das, T47, Apr 2013
%

%if nargin == 0;  help( mfilename ); return; end;

if ~exist('line_pos','var') || isempty( line_pos )
  xlimits = get(gca,'xlim');
  line_pos = round([xlimits(1):xlimits(2)]);
end;
if ~exist('color_code','var'); color_code = 'k'; end;
if ~exist('line_width','var'); line_width = 0.5; end;
if ~exist('low_bd','var'); low_bd = 1; end;
if ~exist('high_bd','var'); high_bd = 1; end;

ylim = get(gca,'ylim');
if low_bd == 1; y_bd(1) = ylim(1); end;
if low_bd == 0; y_bd(1) = 0; end;
if high_bd == 1; y_bd(2) = ylim(2); end;
if high_bd == 0; y_bd(2) = 0; end;
  
hold on;
for i = 1:length( line_pos );
  hold on;
  plot( 0.5+line_pos(i)*[1 1], [y_bd(1) y_bd(2)],'-','color',color_code,'linewidth',line_width); 
end;
hold off;
