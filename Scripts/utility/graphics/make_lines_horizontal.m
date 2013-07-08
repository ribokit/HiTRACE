function make_lines_horizontal(line_pos, colorcode, linewidth, linestyle)
%
% make_lines_horizontal(line_pos,colorcode,linewidth);
%

%if nargin == 0;  help( mfilename ); return; end;

if ~exist( 'line_pos', 'var' ) || isempty( line_pos );
    ylimits = get(gca,'ylim');
    line_pos = round( [ylimits(1):ylimits(2)] );
end;

if ~exist('colorcode', 'var'); colorcode = 'k'; end;
if ~exist('linewidth', 'var'); linewidth = 1; end;
if ~exist('linestyle', 'var'); linestyle='-'; end;

xlim = get(gca, 'xlim');
hold on;
for i = 1:length( line_pos );
  hold on;
  plot( [xlim(1) xlim(2)], 0.5+line_pos(i)*[1 1], '-',...
	'LineWidth', linewidth,...       
	'Color', colorcode,...
    'LineStyle',linestyle );
end;
hold off;
