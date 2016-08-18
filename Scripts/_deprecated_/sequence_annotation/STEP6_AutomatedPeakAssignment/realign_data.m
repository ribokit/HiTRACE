function d_align = realign_data( d, xsel_all,refcol);

if ~exist('refcol')
  refcol = 1;
end

xsel = xsel_all(:,refcol);

num_lanes = size( d, 2 );
num_pixels = size( d, 1 );

for n = 1:num_lanes
  new_x = interp1( xsel_all(:,refcol), xsel_all(:,n), [1:num_pixels],'linear','extrap');
  d_align(:,n) = interp1( 1:num_pixels, d(:,n), new_x,'linear','extrap' );
end
image( d_align )

ymin = round(min(min( xsel))) - 50;
ymax = round(max(max( xsel))) + 50;

ylim( [ ymin ymax ] );

%set(gca,'ytick',xsel,'yticklabel',[] );
hold on
for k = 1:length( xsel )
  plot( [0.5 num_lanes+0.5], xsel(k)*[1 1],'r' );
end
hold off
