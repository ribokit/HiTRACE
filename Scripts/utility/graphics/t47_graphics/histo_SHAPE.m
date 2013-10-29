function histo_SHAPE(d, name)

if ~exist('name','var'); name = ''; end;

nbins = floor(length(d)/5);
if length(d) > 100;
    nbins = nbins*2;
end;

[count,bins] = hist(d,nbins);

figure(); clf;
bar(bins,count/length(d)*100,'r');
legend(name);
ylabel('Count Percentage');
xlabel('SHAPE reactivity');

fprintf('"mean" = %d.\n',mean(d(d<mean(d))));
fprintf('median = %d.\n',median(d));
