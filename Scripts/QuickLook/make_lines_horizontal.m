function make_lines_horizontal(line_pos,colorcode,linewidth);

if ~exist('colorcode')
  colorcode='r';
end
if ~exist('linewidth')
  linewidth = 1;
end

xlim = get(gca,'xlim');
hold on
for i = 1:length( line_pos );
  hold on
  plot( [xlim(1) xlim(2)], 0.5+line_pos(i)*[1 1], '-',...
	'linewidth',linewidth,...
	'color',colorcode); 
end
hold off
