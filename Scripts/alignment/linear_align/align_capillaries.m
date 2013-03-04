function data_align = align_capillaries( data , refcol, reflane, REFINE);
% ALIGN_CAPILLARIES:  Data from ABI readin -- align based on channel 'refcol'.
%
%   data_align = align_capillaries( data , refcol, reflane );
%
% (C) R. Das & S.R. Yoon, 2009-2011
%

if nargin == 0;  help( mfilename ); return; end;

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
  % no longer in use!
  %[d_align, x_realign] = align_to_first_REFINE( d, 0, reflane );
  fprintf( 'ERROR! REFINE is no longer supported!\n');
  return;
else
  % what is the difference between ver2 and ver3?
  [d_align, x_realign] = align_to_first_ver2( d, 0, reflane );
end

x = [1:numpts]';

data_shift=zeros(numpts, 4); % sryoon
data_align{num_capillaries} = data_shift; % sryoon

for i = 1:num_capillaries
  data_shift = interp1( [1:length(data{i})]', data{i}, x_realign(:,i), 'linear',0.0 ); % sryoon
  data_align{i} = data_shift;
end


