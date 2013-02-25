function area_bsub = background_subtracter_LP( area_peak, backgd_cols, normbins, area_pred );

area_peak = quick_norm( area_peak, normbins );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(7)
clf;

goodbins = normbins;
backgd_estimate = mean( area_peak(:, backgd_cols), 2 );

which_col_to_fit = [1:size(area_peak,2)];
area_bsub =[];

for  i = which_col_to_fit;

  s = area_peak(goodbins,i);
  y = area_pred(goodbins,i);
  b = backgd_estimate( goodbins );

  params = background_fit_LP( s, y, b );

  fprintf( '%3d. Background norm:  %8.3f.  Signal strength: %8.3f\n', i, params(1),params(2) );

  area_bsub(:,i) = area_peak(:,i) - params(1)*backgd_estimate;
  %if  ( params(2) > 0.0 )
  %  area_bsub(:,i) = area_bsub(:,i)/params(2);
  %end
  if  ( params(1) > 0.0 )
    area_bsub(:,i) = area_bsub(:,i)/params(1);
  end
  

  %pause;
end

clf;
subplot(1,3,1)
image( 20*area_peak );

subplot(1,3,2)
image( 20 * area_bsub );

subplot(1,3,3)
image( 100*area_pred );

colormap( 1- gray(100));


