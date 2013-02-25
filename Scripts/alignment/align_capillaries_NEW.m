function data_align = align_capillaries_NEW( data , refcol, reflane);

if ~exist( 'refcol')
  refcol = 4;
end
if ~exist( 'reflane')
  reflane = 1;
end

num_capillaries = size( data, 2 );
numpts = size( data{1}, 1);

d = zeros(numpts, num_capillaries); % sryoon
for i = 1:num_capillaries
  d(:,i) = data{i}(:,refcol);
end

[d_align, x_realign] = align_to_first_NEW( d, 0, reflane );

x = [1:numpts]';

data_shift=zeros(numpts, 4); % sryoon
data_align{num_capillaries} = data_shift; % sryoon
for i = 1:num_capillaries
%   for m = 1:4
%     d  = data{i}(:,m);
%     data_shift(:,m) = interp1( x, d, x_realign(:,i),'linear',0.0 );
%   end
  data_shift = interp1(x, data{i}, x_realign(:,i), 'linear',0.0 ); % sryoon
  data_align{i} = data_shift;
end

