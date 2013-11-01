function [data_align, x_realign] = align_to_first_ver2( data, PLOT_STUFF, refcol );
%
% ALIGN_TO_FIRST_VER2:  (linear-time) alignment of a matrix of electropheretic traces to first trace
%
%   [data_align, x_realign] = align_to_first_ver2( data, PLOT_STUFF, refcol )
%
%  Optimizes correlation coefficient by grid search + Fast Fourier Transform
%
%  Uses 'peakified' traces to help prevent outliers from distorting alignment.
%
% (C) R. Das & S.R. Yoon, 2009-2011
%

if nargin == 0;  help( mfilename ); return; end;

FULL_SIGNAL_WINDOW = 5000;

if ~exist('refcol', 'var');  refcol = 1; end

num_capillaries = size( data,2 );

if ~exist('PLOT_STUFF', 'var');  PLOT_STUFF = 0; end

% main difference with align_to_first_ver3 -- use of 'peakify'.
peak = peakify( data );
data = peak; % jk

%d_ref = baseline_subtract( data(:,refcol) );
d_ref = data(:,refcol);

[minbin_ref, middlebin_ref, maxbin_ref] = get_signal_bins( d_ref, FULL_SIGNAL_WINDOW );
d1 = extract_profile( d_ref , minbin_ref, maxbin_ref );

numpts_d_ref = length( d_ref ); % sryoon
x = [ 1: numpts_d_ref ]; % sryoon
x_realign = zeros(numpts_d_ref,num_capillaries); % sryoon
data_align = zeros(numpts_d_ref,num_capillaries); % sryoon

numpts = length(d1); % sryoon
%shifts = [(-numpts+1):(numpts-1)]; % sryoon
max_shift = 200;
shifts = [-max_shift:max_shift]; % rhiju
scales = [0.95:0.005:1.05]; % sryoon

if parallelization_exists()
    parfor n = 1: num_capillaries
        fprintf(1,'Calculating correlation...%d\n',n);
        [data_realign(:,n), x_realign(:,n) ] = align_to_first_inner_loop( data(:,n), FULL_SIGNAL_WINDOW, scales, shifts, d1, x, minbin_ref, maxbin_ref, PLOT_STUFF );
    end
else
    fprintf('\n'); revStr = ' '; fprintf(' \n');
    
    for n = 1: num_capillaries
        revStr = lprintf(revStr,['Calculating correlation ', num2str(n), ' of ', num2str(num_capillaries), ' ... \n'], 2);
        [data_realign(:,n), x_realign(:,n) ] = align_to_first_inner_loop( data(:,n), FULL_SIGNAL_WINDOW, scales, shifts, d1, x, minbin_ref, maxbin_ref, PLOT_STUFF );
    end
end


data_align = data_realign;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data_align, x_realign] =  align_to_first_inner_loop( d, FULL_SIGNAL_WINDOW, scales, shifts, d1, x,minbin_ref, maxbin_ref, PLOT_STUFF );

%d     = baseline_subtract( data(:,n) );

[minbin, middlebin, maxbin ] = get_signal_bins( d, FULL_SIGNAL_WINDOW );
%d2 = extract_profile( dezinger(d), minbin, maxbin );

d2 = extract_profile( d, minbin, maxbin );
%d2 = d2 - smooth( d2, 100);

%numpts = length(d1); % sryoon

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%shifts = [(-numpts+1):(numpts-1)]; %sryoon

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
    
    scales = best_scale + [-0.01:0.0005:0.01];
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
x_realign = new_x;
%da = baseline_subtract( da );
data_align = da;

if PLOT_STUFF
    plot( x, d_ref/mean(d_ref(minbin_ref:maxbin_ref)), 'b', ...
        x, da/mean(da( minbin_ref:maxbin_ref)), 'r' );
    axis([ minbin-200 maxbin+200 -10 40]);
    pause;
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  [minbin, middlebin, maxbin ] = get_signal_bins( d, FULL_SIGNAL_WINDOW );
%Find continuous chunk with the most signal.

% d = baseline_subtract( d );
% d = d - smooth(d, 100 );

%numpts = length( d ); %sryoon

% d_sum = cumsum( abs(d) );
% d_sum2 = circshift( d_sum, FULL_SIGNAL_WINDOW );
% tot_signal = d_sum - d_sum2;
tot_signal = d;

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

%minbin = 1500; %jk
%maxbin = 4000; %jk

minbin = 101;
maxbin = length( d )-100;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d1 = extract_profile( d, minbin, maxbin );

d1 = d( minbin:maxbin) - mean( d(maxbin+[1:100]) );

% CUTOFF = 1000;
CUTOFF = 200000 * abs(mean(d1)) / (1+std(d1)^2);
d1 = ( min(max( d1, 0 ),CUTOFF) ).^0.5;

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
best_scale = scales( best_scale_index );
best_shift = shifts( best_shift_index );
return;