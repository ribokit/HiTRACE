function [bandpos, band] = xsel_label_split(xsel, seqpos, sequence, offset, H, h_length, h_sp)

% read in band names (Y-axis)
bandpos = cell(length(seqpos), 2);
for i = 1:length(seqpos)
    bandpos{i, 1} = [sequence(seqpos(i) - offset), num2str(seqpos(i))];
    bandpos{i, 2} = xsel(i);
end;

% split band position array into h*h_length array
band = cell(H, size(bandpos, 1), 2);
for i = 1:H
    ymin = h_length * (i - 1) + 1; ymax = h_length * i;
    for j = 1:size(bandpos, 1)
        if bandpos{j, 2}>=ymin && bandpos{j, 2}<=ymax;
            band{i, j, 1} = bandpos{j, 1}; 
            band{i, j, 2} = bandpos{j, 2} + h_sp - (ymin - 1);
        end;
    end;
end;
