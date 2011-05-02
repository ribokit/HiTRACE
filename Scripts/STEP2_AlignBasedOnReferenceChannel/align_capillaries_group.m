function data_align = align_capillaries_group( data, refcol ,reflane, refgrp)
% by kprotoss, group align

if ~exist( 'refcol')
  refcol = 4;
end
if ~exist( 'reflane')
  reflane = 1;
end
if ~exist('refgrp')
    refgrp = 1;
end

num_grp = length( data );
numpts = size( data{1}{1}, 1);
num_capillaries = 0;

%find group x_align
d = zeros(numpts, num_grp); 
for i = 1:num_grp
  v = [1 round((1+size(data{i},2))/2) size(data{i},2)];
  d(:,i) = data{i}{v(1)}(:,refcol)+data{i}{v(2)}(:,refcol)+data{i}{v(3)}(:,refcol);
  num_capillaries = num_capillaries + length(data{i});
end

%[d_align, x_realign] = align_to_first_OLD( d, 0, refgrp);
%[d_align, x_realign] = align_to_first_ver2( d, 0, refgrp);
[d_align, x_realign] = align_to_first_ver3( d, 0, refgrp);

% re-align to first group.
x = [1:numpts]';

data_shift=zeros(numpts, 4);
data_align{num_capillaries} = data_shift;
count =0;
for i = 1:num_grp
    for j = 1:length(data{i})
        count = count + 1;
        data_shift = interp1(x, data{i}{j}, x_realign(:,i), 'linear',0.0 );
        data_align{count} = data_shift;
    end
end