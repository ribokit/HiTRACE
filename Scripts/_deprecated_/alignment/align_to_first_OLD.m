function [data_align, x_realign] = align_to_first_OLD( data, PLOT_STUFF, refcol );

FULL_SIGNAL_WINDOW = 2000;

if ~exist('refcol')
  refcol = 1;
end

num_capillaries = size( data,2 );
%numpts = size( data, 1); % sryoon

if ~exist('PLOT_STUFF')
  PLOT_STUFF = 0;
end

d_ref = baseline_subtract( data(:,refcol) );
[minbin_ref, middlebin_ref, maxbin_ref] = get_signal_bins( d_ref, FULL_SIGNAL_WINDOW );
%d1 = extract_profile( dezinger(d_ref) , minbin_ref, maxbin_ref );
d1 = extract_profile( d_ref , minbin_ref, maxbin_ref );
%d1 = d1 - smooth( d1, 100);

%BEGIN sryoon
numpts_d_ref = length( d_ref );
x = [ 1: numpts_d_ref ]; 
x_realign = zeros(numpts_d_ref,num_capillaries); 
data_align = zeros(numpts_d_ref,num_capillaries);
numpts = length(d1); 
%shifts = [(-numpts+1):(numpts-1)]; 
max_shift = 500;
shifts = [-max_shift:max_shift]; % rhiju
scales = [0.93:0.005:1.08]; 
%END sryoon

for n = 1: num_capillaries
  d     = baseline_subtract( data(:,n) );

  [minbin, middlebin, maxbin ] = get_signal_bins( d, FULL_SIGNAL_WINDOW );
  %d2 = extract_profile( dezinger(d), minbin, maxbin );
  d2 = extract_profile( d, minbin, maxbin );
  %d2 = d2 - smooth( d2, 100);
  
  %numpts = length(d1); % sryoon

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %shifts = [(-numpts+1):(numpts-1)]; %sryoon

  fprintf(1,'Calculating correlation...%d\n',n);

  % First do grid search
  %scales = [0.95:0.005:1.05]; %sryoon
  %scales = [0.8:0.01:1.2];
  [best_scale, best_shift] = find_best_scale_shift( scales, shifts, ...
						    d1, d2 );

  % Following helps a little, but doesn't really allow us to skip a
  % fine grid search!
  REFINE = 0;
  if REFINE
    scales = best_scale + [-0.02:0.005:0.02];
    [best_scale, best_shift] = find_best_scale_shift( scales, shifts, ...
						      d1, d2 );

    scales = best_scale + [-0.01:0.001:0.01];
    [best_scale, best_shift] = find_best_scale_shift( scales, shifts, ...
						      d1, d2 );
  end
  
  %  Then refine
  %  [Doesn't help too much.]
  MINIMIZE_REFINE = 0;
  if MINIMIZE_REFINE
    options = optimset( 'MaxFunEvals',10,'Display','off','TolX',0.001);
    best_scale = ...
	fminbnd( 'best_align', best_scale - 0.02, best_scale + 0.02, options, d1, d2 );
    [dummy, best_shift_index] = best_align( best_scale, d1, d2 );
    best_shift = shifts( best_shift_index );
  end
  
  if PLOT_STUFF
    subplot(1,1,1);    
    if 0
      d2a = interp1( x, d2, scales(best_scale_index)*x, 'linear',0);
      %plot( x, (d1.^2)/ mean( d1.^2), 'b',x, (d2a.^2)/mean(d2.^2),'r' );
      plot( x, d1, 'b',x+shifts( best_shift_index), d2a,'r' );
      axis([0 numpts 0 20]);
      pause;    
    end
  end
  
  %numpts = length( d_ref ); % sryoon
  %x = [ 1: numpts_d_ref ]; % sryoon
  new_x = best_scale * (x-minbin_ref+1 - best_shift)  + minbin - 1  ;
  da = interp1( x, d, new_x, 'linear',0);
  x_realign(:,n) = new_x;
  da = baseline_subtract( da );
  data_align(:,n) = da;

  if PLOT_STUFF
    plot( x, d_ref/mean(d_ref(minbin_ref:maxbin_ref)), 'b', ...
	  x, da/mean(da( minbin_ref:maxbin_ref)), 'r' );
    axis([ minbin-200 maxbin+200 -10 40]);    
    pause;
  end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  [minbin, middlebin, maxbin ] =get_signal_bins( d, FULL_SIGNAL_WINDOW );
%Find continuous chunk with the most signal.

d = baseline_subtract( d );
d = d - smooth(d, 100 );

%numpts = length( d ); %sryoon

d_sum = cumsum( abs(d) );
d_sum2 = circshift( d_sum, FULL_SIGNAL_WINDOW );
tot_signal = d_sum - d_sum2;

[dummy, maxbin] = max( tot_signal((1+FULL_SIGNAL_WINDOW):end) );
maxbin = maxbin + FULL_SIGNAL_WINDOW;
minbin = maxbin - FULL_SIGNAL_WINDOW + 1;
middlebin = floor( 0.75*maxbin + 0.25*minbin  );

if 0 
clf
plot( d*100,'r' );
hold on
plot( tot_signal );
hold off
pause;
end

minbin = 1500;
maxbin = 4500;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d1 = extract_profile( d, minbin, maxbin );

d1 = d( minbin:maxbin) - mean( d(maxbin+[1:100]) );

CUTOFF = 500;
d1 = ( min(max( abs(d1), 0 ),CUTOFF) ).^0.5;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function   [best_scale, best_shift] = find_best_scale_shift( scales, ...
						  shifts, d1, d2 );

% numscales = length( scales );
% x = 1:length(d1); 
% d2x = [];
%  for k = 1:numscales
%    d2x(:,k) = interp1( x, d2, scales(k)*x,'linear', d2(end) );
%  end
% conv_matrix = conv2( d1, d2x(end:-1:1,:));

% BEGIN sryoon
len_d1 = length(d1);
n_fft_row = 2*len_d1 - 1;
n_fft_col = length(scales);
x = (1:len_d1)';
d2x=interp1(x, d2, x*scales, 'linear',d2(end));
conv_matrix = real(ifft2(fftn(d1,[n_fft_row n_fft_col]) ...
    .* fftn(d2x(end:-1:1,:), [n_fft_row n_fft_col])));
% END sryoon

shifts_all = [ -(len_d1-1):(len_d1-1) ];
goodpoints = find( shifts_all >= min(shifts) & shifts_all <= max(shifts));
conv_matrix = conv_matrix(goodpoints,:);

[ dummy, best_scale_index ] = max( max( conv_matrix ) );
[ dummy, best_shift_index ] = max( conv_matrix(:,best_scale_index) );
best_scale = scales( best_scale_index);
best_shift = shifts( best_shift_index );

