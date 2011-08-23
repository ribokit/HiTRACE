function make_logL_contour_plot( logLout, param1, param2, param1_name, param2_name )

if ~exist( 'param1_name' ); param1_name = 'K'; end;
if ~exist( 'param2_name' ); param2_name = 'n_{Hill}'; end;
  
contours = max(max(logLout)) - [0:30, 30:10:200];


[x_axis, x_is_linear, min_x, max_x, xtick, x_name] = get_lin_or_log_axis( param1, param1_name );
[y_axis, y_is_linear, min_y, max_y, ytick, y_name] = get_lin_or_log_axis( param2, param2_name );

contour( x_axis, y_axis, -logLout',-contours);

xlim( [ min_x max_x ] );
ylim( [ min_y max_y ] );

set(gca,'xgrid','on','ygrid','on', 'xtick', xtick, 'ytick',ytick);
xlabel( x_name );
ylabel( y_name );

%%%%%%%%%%%%%%%%%%%%%
function [x_axis, x_is_linear, min_x, max_x, xtick, x_name] = get_lin_or_log_axis( param, param_name );

x_is_linear = 0.0;
  
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
else
  x_axis = log( param )/ log(10);
  min_x = min( x_axis );
  max_x = max( x_axis );
  xtick = [floor(min_x) : 0.5:ceil( max_x ) ];
  x_name = ['log_{10} ',  param_name ];
end
