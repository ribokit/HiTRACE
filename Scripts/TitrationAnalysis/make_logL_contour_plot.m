function make_logL_contour_plot( logLout, param1, param2, param1_name, param2_name, p1_best, p2_best )

if ~exist( 'param1_name' ); param1_name = 'K'; end;
if ~exist( 'param2_name' ); param2_name = 'n_{Hill}'; end;
  
cont_levels = [0:2:20, 25, 30, 35, 40:20:200];
contours = max(max(logLout)) - cont_levels;
fprintf( ['Making contours lower than max log-likelihood point by: ', num2str(cont_levels ),'\n'] );

if ~exist( 'p1_best' ) | ~exist( 'p2_best' ); p1_best = []; p2_best = []; end;

[x_axis, x_is_linear, min_x, max_x, xtick, x_name, x_best] = get_lin_or_log_axis( param1, param1_name, p1_best );
[y_axis, y_is_linear, min_y, max_y, ytick, y_name, y_best] = get_lin_or_log_axis( param2, param2_name, p2_best );

contour( x_axis, y_axis, -logLout',-contours);

if ~isempty( x_best) & ~isempty( y_best )
  hold on; plot( x_best, y_best, 'ko','markerfacecolor','k','linewidth',3,'markersize',6 ); hold off;
end

xlim( [ min_x max_x ] );
ylim( [ min_y max_y ] );

set(gca,'xgrid','on','ygrid','on', 'xtick', xtick, 'ytick',ytick,'linew',2,'fontsize',14,'fontw','bold');
xlabel( x_name );
ylabel( y_name );

%%%%%%%%%%%%%%%%%%%%%
function [x_axis, x_is_linear, min_x, max_x, xtick, x_name, x_best] = get_lin_or_log_axis( param, param_name, p_best );

x_is_linear = 0.0;
x_best = [];
  
if (length( param ) < 3  | ...
    abs( (param(3) - param(2)) / (param(2)-param(1)) - 1 ) < 0.001 )
  x_is_linear = 1;
end

if x_is_linear
  x_axis = param;
  min_x = 0;
  max_x = max( x_axis );
  xtick = [0:0.5:ceil( max_x )];
  x_name = param_name;
  if ~isempty( p_best ); x_best = p_best; end;
else
  x_axis = log( param )/ log(10);
  min_x = min( x_axis );
  max_x = max( x_axis );
  xtick = [floor(min_x) : 0.5:ceil( max_x ) ];
  x_name = ['log_{10} ',  param_name ];
  if ~isempty( p_best ); x_best = log( p_best )/log(10); end;
end

