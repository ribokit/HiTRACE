function make_image( d, d_unmod )
subplot(1,2,1); 
image( d );
title( 'initial');

subplot(1,2,2); 
image( d_unmod );
title( 'after correction');

colormap( 1- gray(100 ) );
