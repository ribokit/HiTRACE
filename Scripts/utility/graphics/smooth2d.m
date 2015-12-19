function x_smooth = smooth2d( x, niter );
if ~exist( 'niter' )  niter = 2; end
B = [ 0 0.1 0; 0.1 0.6 0.1; 0 0.1 0];
x_smooth = conv2( x, B, 'same' );

for  n = 1:niter-1;
    x_smooth = conv2( x_smooth, B, 'same' );
end