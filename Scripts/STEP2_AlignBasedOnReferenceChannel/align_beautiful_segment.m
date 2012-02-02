function [d_realign,da_realign] =  align_beautiful_segment( d, da, align_pts, refcol );

if ~exist( 'align_pts' )
  align_pts = [1:size(d,1)];
end
if ~exist('refcol')
  refcol = 1;
end

d_to_realign  = d( align_pts,: );
da_to_realign = da( align_pts,: );

numpts = size( d_to_realign, 1);
num_capillaries = size( d_to_realign, 2 );
[d_align,x_realign] = align_to_first_ver2( 10*da_to_realign, 0, refcol );

x = [1:numpts];
d_realign = [];
for i = 1:num_capillaries
  d_realign(:,i) = interp1(x, d_to_realign(:,i), x_realign(:,i), 'linear',0.0 );
  da_realign(:,i) = interp1(x, da_to_realign(:,i), x_realign(:,i), 'linear',0.0 );
end

subplot(1,2,1);
image( 40*d_to_realign );
subplot(1,2,2);
image( 40*d_realign );
colormap( 1- gray(100) );