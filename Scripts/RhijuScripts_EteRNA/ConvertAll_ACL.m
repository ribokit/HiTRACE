PREFIX = 'NewRDAT';
filelist = ls('OldRDAT\*.rdat');

for i = 1:size(filelist,1)
    convert_024_030(['OldRDAT\', filelist(i,:)],['NewRDAT\', filelist(i,:)]);
end