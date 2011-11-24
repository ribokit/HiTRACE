% by kprotoss
% find time range automatically

function [ymin ymax] = findTimeRange(d)
for j = 1:size(d,2);
    y = edge(d(:,j));
    y = find( y ~= 0);

    min_list(j) = y(1);
    max_list(j) = y(end);
end

ymin = round(median(min_list) - 100);
ymax = round(median(max_list) + 100);