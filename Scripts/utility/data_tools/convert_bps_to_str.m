function bps_str = convert_bps_to_str(bps)

bps_str = cell(size(bps,1),1);
for i = 1:size(bps,1);
    bps_str{i} = [num2str(bps(i,1)),';',num2str(bps(i,2))];
end;