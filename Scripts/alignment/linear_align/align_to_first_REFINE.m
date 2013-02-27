function [data_align, x_realign] = align_to_first_REFINE( data, PLOT_STUFF, refcol );

if ~exist('refcol')
  refcol = 1;
end

num_capillaries = size( data,2 );

if ~exist('PLOT_STUFF')
  PLOT_STUFF = 0;
end


minbin = 1500;
maxbin = 8000;
minbin_ref = 1500;
maxbin_ref = 8000;

d_ref = baseline_subtract( data(:,refcol) );
d1 = extract_profile( d_ref , minbin_ref, maxbin_ref );

%BEGIN sryoon
numpts_d_ref = length( d_ref );
x = [ 1: numpts_d_ref ]; 
x_realign = zeros(numpts_d_ref,num_capillaries); 
data_align = zeros(numpts_d_ref,num_capillaries);
numpts = length(d1); 
shifts = [(-numpts+1):(numpts-1)]; 
scales = [0.95:0.005:1.05]; 
%END sryoon
quad_coeffs = [ -1e-5: 1e-6: 1e-5];

NUM_WINDOWS = 10;

for n = [1: num_capillaries]

  d     = baseline_subtract( data(:,n) );
  d2 = extract_profile( d, minbin, maxbin );

  fprintf(1,'Calculating correlation...%d\n',n);

  % First do grid search -- linear
  [best_scale, best_shift] = find_best_scale_shift( scales, shifts, ...
						    d1, d2 );    
  
  % Next do search for quadratic.
  [best_quad_coeff, best_shift] = find_best_scale_shift_quadratic( best_scale, quad_coeffs, shifts, ...
						    d1, d2 );    
  
  x_shift = x - minbin_ref + 1 - best_shift;
  new_x = best_scale * x_shift  +  x_shift.*x_shift*best_quad_coeff   + minbin - 1  ;
  da = interp1( x, d, new_x, 'linear',0);
  x_realign(:,n) = new_x;
  da = baseline_subtract( da );
  data_align(:,n) = da;

  
  if PLOT_STUFF
    plot( x, d_ref/mean(d_ref(minbin_ref:maxbin_ref)), 'b', ...
	  x, da/mean(da( minbin_ref:maxbin_ref)), 'r' );
    axis([ minbin-200 maxbin+200 -10 40]);    
  end

end


return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d1 = extract_profile( d, minbin, maxbin );

d1 = d( minbin:maxbin) - mean( d(maxbin+[1:100]) );

CUTOFF = 500;
d1 = ( min(max( d1, 0 ),CUTOFF) ).^0.5;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function   [best_scale, best_shift] = find_best_scale_shift( scales, ...
						  shifts, d1, d2 );
len_d1 = length(d1);
n_fft_row = 2*len_d1 - 1;
n_fft_col = length(scales);
x = (1:len_d1)';
d2x=interp1(x, d2, x*scales, 'linear',d2(end));
conv_matrix = real(ifft2(fftn(d1,[n_fft_row n_fft_col]) ...
    .* fftn(d2x(end:-1:1,:), [n_fft_row n_fft_col])));

[ dummy, best_scale_index ] = max( max( conv_matrix ) );
[ dummy, best_shift_index ] = max( conv_matrix(:,best_scale_index) );
best_scale = scales( best_scale_index);
best_shift = shifts( best_shift_index );

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [best_quad_coeff, best_shift] = find_best_scale_shift_quadratic( ...
    scale, quad_coeffs,...
    shifts, d1, d2 );

len_d1 = length(d1);
n_fft_row = 2*len_d1 - 1;
n_fft_col = length(quad_coeffs);
x = (1:len_d1)';
d2x=interp1(x, d2, x * ( scale*ones(1,n_fft_col)) + ...
	    (x.*x)*quad_coeffs, 'linear',d2(end));

%plot( d2x );
%hold on
%plot( d1,'linew',2);
%hold off
%pause;

conv_matrix = real(ifft2(fftn(d1,[n_fft_row n_fft_col]) ...
    .* fftn(d2x(end:-1:1,:), [n_fft_row n_fft_col])));


%image( conv_matrix*0.0001 )
%pause
[ dummy, best_quad_index ] = max( max( conv_matrix ) );
[ dummy, best_shift_index ] = max( conv_matrix(:,best_quad_index) );
best_quad_coeff = quad_coeffs( best_quad_index);
best_shift = shifts( best_shift_index );

