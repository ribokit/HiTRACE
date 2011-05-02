function data_align = align_capillaries_ver2( data , refcol, reflane, ref );

if ~exist( 'refcol')
  refcol = 4;
end
if ~exist( 'reflane')
  reflane = 1;
end

%jk begin
data_orig = data;
data = data_cleaning(data);
%jk end

num_capillaries = size( data, 2 );
numpts = size( data{1}, 1);
%colorcode = jet( num_capillaries ); % sryoon

d = zeros(numpts, num_capillaries); % sryoon
for i = 1:num_capillaries
  d(:,i) = data{i}(:,refcol);
end

[d_align, x_realign] = align_to_first_ver2_OLD( d, 0, reflane );

x = [1:numpts]';

data_shift=zeros(numpts, 4); % sryoon
data_align{num_capillaries} = data_shift; % sryoon
for i = 1:num_capillaries
%   for m = 1:4
%     d  = data{i}(:,m);
%     data_shift(:,m) = interp1( x, d, x_realign(:,i),'linear',0.0 );
%   end
  data_shift = interp1(x, data{i}, x_realign(:,i), 'linear',0.0 ); % sryoon %jk
  data_align{i} = data_shift;
end


% Scaling
if 0
d = zeros(numpts, num_capillaries); % sryoon
for i = 1:num_capillaries
  d(:,i) = data_align{i}(:,refcol);
  d(2000:2500,i) = d(2000:2500,i).*100;
  d(1:2000,i) = 0;%d(1:2000,i)./2;
  d(2500:end,i) = 0;%d(2500:end,i)./2;
end

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
end