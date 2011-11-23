function data_align = align_capillaries( data , refcol, reflane, REFINE);
% ALIGN_CAPILLARIES:  Data from ABI readin -- align based on channel 'refcol'.
%
%   data_align = align_capillaries( data , refcol, reflane, REFINE);
%
% (C) R. Das & S.R. Yoon, 2009-2011
%

if ~exist( 'refcol')
  refcol = 4;
end
if ~exist( 'reflane')
  reflane = 1;
end
if ~exist( 'REFINE')
  REFINE = 0;
end

num_capillaries = size( data, 2 );
numpts = size( data{1}, 1);
%colorcode = jet( num_capillaries ); % sryoon

for i = 1:num_capillaries
  numpts = max( numpts,  length(data{i}(:,refcol) ));
end
d = zeros(numpts, num_capillaries); 
for i = 1:num_capillaries
  numpts = length(data{i}(:,refcol) );
  d(1:numpts,i) = data{i}(:,refcol);
end

if REFINE
  [d_align, x_realign] = align_to_first_REFINE( d, 0, reflane );
else
  %[d_align, x_realign] = align_to_first_OLD( d, 0, reflane );
  [d_align, x_realign] = align_to_first_ver2( d, 0, reflane );
  %[d_align, x_realign] = align_to_first_ver3( d, 0, reflane );
end

x = [1:numpts]';

data_shift=zeros(numpts, 4); % sryoon
data_align{num_capillaries} = data_shift; % sryoon
for i = 1:num_capillaries
%   for m = 1:4
%     d  = data{i}(:,m);
%     data_shift(:,m) = interp1( x, d, x_realign(:,i),'linear',0.0 );
%   end
  data_shift = interp1( [1:length(data{i})]', data{i}, x_realign(:,i), 'linear',0.0 ); % sryoon
  data_align{i} = data_shift;
end


