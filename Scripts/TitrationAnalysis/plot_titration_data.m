function plot_titration_data( data, resnum, conc, ...
		    pred_fit, sigma_at_each_residue, lane_normalization, ...			      
			      conc_fine, pred_fit_fine );

numres  = size( data,1 ); 
numconc = size( data,2 ); 

if exist( 'lane_normalization' )
  for j = 1:length( conc )
    data(:,j) = data(:,j) / lane_normalization(j);
    %pred_fit(:,j) = pred_fit(:,j) / lane_normalization(j);
  end
end

figure(2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clf;
colorcode = jet(numconc);

plot_offset = 0.0;

for j = 1:numconc
  plot( resnum, plot_offset*(j-1) + data(:,j), '.','color',colorcode(j,:),...
	'markerfacecolor',colorcode(j,:));
  hold on
end
plot( resnum, sigma_at_each_residue,'color', [0.5 0.5 0.5] );
for j = 1:numconc
  plot( resnum, plot_offset*(j-1) + pred_fit(:,j), '-','color',colorcode(j,:));
end
xlabel('Residue number');
ylabel('Peak intensity');

for m = 1:length( conc ); legends{m} = num2str( conc(m) ); end;
legends{ length( conc ) + 1 } = 'fitted stdev';
legend( legends );

% outliers?
vals =  max( pred_fit' );
gp = find( vals <= ( mean( vals ) + 3*std( vals ) ) ); % remove outliers
if length(vals(gp) ) > 1
  ylim( [ 0 max( vals( gp ) )]);
else
  ylim( [min( pred_fit )  max(pred_fit)] );
end

xlim( [min(resnum)-1,  max(resnum)+1] );
%set(gca,'ylim',[0 (numconc+1)*plot_offset])
hold off

figure(3)
%subplot(1,2,2);
colorcode = jet(numres);
plot_offset = mean(mean(data));
for i = 1:numres
  semilogx( conc, plot_offset*(i-1) + data(i,:), '.','color',colorcode(i,:),...
	    'markerfacecolor',colorcode(i,:));
  hold on
  plot( conc_fine, plot_offset*(i-1) + pred_fit_fine(i,:), '-','color',colorcode(i,:));

  startpt = max(find(conc>0));
  h = text( conc(startpt), plot_offset*(i-1)+pred_fit(i,startpt), num2str( resnum( i ) ) );
  set(h,'color','k','fontsize',8,'fontweight','bold');
end
ylim2 = [0 (numres+1)*plot_offset+max(max(data))];
if ( size( pred_fit, 1 )== 1  ) ylim2 = [min( data) max(data) ]; end;
set(gca,'ylim',ylim2,'xlim',[ min(conc) max(conc) ])
set(gca,'linew',2,'fontsize',14,'fontw','bold');
xlabel('[M^{2+} (mM)]');
ylabel('Peak intensity');

hold off

