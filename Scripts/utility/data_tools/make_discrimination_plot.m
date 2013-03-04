function make_discrimination_plot( all_TP_lengths, all_FP_lengths, L, input_label );
% make_discrimination_plot( all_TP_lengths, all_FP_lengths, L, input_label );
%
if nargin == 0;  help( mfilename ); return; end;

subplot(2,1,1)
h_TP = hist( all_TP_lengths, L );
h_FP = hist( all_FP_lengths, L );
stairs( L, h_TP, 'b' , 'linewidth', 2 ); hold on
stairs( L, h_FP, 'r' , 'linewidth', 2) ; hold off
xlim( [min(L) max(L) ] );
ylim( [0  max([h_TP h_FP])+1 ]);
set(gca,'linewidth',2);
xlabel( input_label );

accuracy = h_TP./(h_TP+h_FP );

subplot(2,1,2)
plot( L, 100*accuracy, 'o','color',[0.5 0.5 0.5], 'markerfacecolor',[0.5 0.5 0.5] );
errors = sqrt( accuracy.*(1-accuracy)./( h_TP + h_FP ) );
hold on
for i = 1:length(L)
  plot( L(i)* [1 1], 100*( accuracy(i) +  errors(i)* [-1 1] ), 'color', [0.5 0.5 0.5] ); 
end
xlabel( 'Average support value' );
ylabel( 'Accuracy' );

hold on
plot( L, L, 'k' );
set(gca,'linewidth',2);
hold off

axis([0 100 0 100]);
