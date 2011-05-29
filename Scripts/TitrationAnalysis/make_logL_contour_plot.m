function make_logL_contour_plot( logLout, K1, nHill )

contours = max(max(logLout))-[0:30, 30:10:200];
contour( log(K1)/log(10), nHill, -logLout',-contours);
set(gca,'xlim',[-2 2],'ylim',[0 max(nHill)],'xtick',[-2:2],'ytick',[0:0.5:max(nHill)],'xgrid','on','ygrid','on');

ylabel( 'n_{Hill}' );
xlabel( 'log_{10} (K_1/ 1 mM)' );
xlim( log([ min(K1) max(K1)]) /log(10) );
