function make_lines(line_pos,colorcode,linewidth);
%
% make_lines(line_pos,colorcode,linewidth);
%

%if nargin == 0;  help( mfilename ); return; end;

if ~exist('line_pos','var') | isempty( line_pos )
  xlimits = get(gca,'xlim');
  line_pos = round([xlimits(1):xlimits(2)]);
end
if ~exist('colorcode','var')
  colorcode='k';
end
if ~exist('linewidth','var')
  linewidth=1;
end

ylim = get(gca,'ylim');
hold on
for i = 1:length( line_pos );
  hold on
  plot( 0.5+line_pos(i)*[1 1], [ylim(1) ylim(2)],'-','color',colorcode,'linewidth',linewidth); 
end
hold off
