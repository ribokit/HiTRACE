function [ log_L, C_state, input_data_rescale, conc_fine, pred_fit_fine_rescale ] = analyze_titration_with_likelihood( input_data, conc, resnum, param1, param2, whichres, fit_type, C_state_in, plot_res )
%  log_L  = analyze_titration_with_likelihood( input_data, conc, resnum, K1_conc, param2, whichres );
%
% Likelihood-based analysis of structure mapping titration -- optimizes lane normalizaton and calculates
%  errors at each residue while doing a grid search over midpoints and apparent Hill coefficients.
%
%  Inputs:
%  input_data = matrix of input structure mapping data, approximately normalized (e.g., by mean intensity in each lane)
%                 must have dimensions of (number of residues x number of lanes ).
%  conc    = concentration of chemical in titration (e.g., [adenine], or [Mg2+] )
%  resnum  = your favorite numbering of residues
%  K1_conc = concentrations to search over. (Default: 10.^[-3.0 : 0.1 :3.0]).
%  param2   = apparent Hill coefficients to search over. (Default: param2 = [0.0:0.05:4] )
%  whichres = [optional!] just look at this subset of input residues.
%  fit_type = string specifying functional form to search: "hill", "double_exp", "one_two"
%
% (C) Rhiju Das, 2008-2011
%

if ~exist( 'param1' ) | isempty( param1 ); param1 = 10.^[-3.0 : 0.1 :3.0]; end
if ~exist( 'param2' ) | isempty( param2 ); param2 = [0.0:0.05:4]; end;
if ~exist( 'resnum' ) | isempty( resnum ); resnum = [1:size( input_data, 1 ) ]; end;
if exist( 'whichres' ) & ~isempty( whichres )
  for k = 1:length(whichres)
    res_to_fit(k) = find( resnum == whichres(k) );
  end
  input_data = input_data( res_to_fit, : );
  resnum = resnum( res_to_fit );
end
if ~exist( 'fit_type' ); fit_type = 'hill'; end;
if ~exist( 'C_state_in' ); C_state_in = []; end;

log_L_best = -99999999999;
p1_best = 999;
p2_best = 999;

[ f, p1_name, p2_name ] = feval( fit_type, conc, param1(1), param2(1) );

for i = 1: length( param1)
  for j = 1: length( param2) % this could be parallelized for speed..

    p1 = param1( i );
    p2 = param2( j );
    
    PLOTSTUFF = 0;
    [logL, pred_fit, lane_normalization, sigma_at_each_residue] = ...
	do_new_likelihood_fit( input_data, conc, p1, p2, fit_type, [], C_state_in );
    
    log_L(i,j) = logL;
    sigma_all(:,i,j) = [sigma_at_each_residue' std( lane_normalization )];
        
    fprintf(1,'%s %8.4f. %s %8.4f ==>  LogL %8.4f\n', p1_name, p1, p2_name, p2, log_L(i,j));
    
    if (log_L(i,j) > log_L_best)
      log_L_best = log_L(i,j);
      p1_best = p1;
      p2_best = p2;
    end

%    sigma_at_each_residue'
%    plot_data( input_data, resnum, conc, pred_fit, sigma_at_each_residue );
%    pause;

  end;
end

[logLbest, pred_fit, lane_normalization, sigma_at_each_residue, C_state ] =...
    do_new_likelihood_fit( input_data, conc, p1_best, p2_best, fit_type, [],  C_state_in );
%lane_normalization

% let's try to do a good job of figuring out errors by looking over likelihood
hold on
if ( length( param1 ) > 1 & length( param2 ) > 1 )
 c = contour( log_L, log_L_best - 2.0 );
 p1_low    = interp1( 1:length(param1), param1,  min( c(2,2:end) ) );
 p1_high   = interp1( 1:length(param1), param1,  max( c(2,2:end) ) );
 p2_low    = interp1( 1:length(param2), param2,  min( c(1,2:end) ) );
 p2_high   = interp1( 1:length(param2), param2,  max( c(1,2:end) ) );
else
 [p1_low,p1_high] = get_1D_error( param1, log_L );
 [p2_low,p2_high] = get_1D_error( param2, log_L );
end
titlestring = sprintf( '%s = %6.4f + %6.4f - %6.4f\n', p1_name, p1_best, p1_high-p1_best, p1_best - p1_low );
titlestring = [titlestring, sprintf( '%s = %4.2f + %4.2f - %4.2f', p2_name, p2_best, p2_high - p2_best, p2_best - p2_low ) ];
fprintf( [titlestring,'\n'] )


log_10_conc = log( conc( find( conc > 0 ) ) ) / log( 10 );
conc_fine = 10.^[min(log_10_conc):0.01:max(log_10_conc)];
f = feval( fit_type, conc_fine, p1_best, p2_best);
pred_fit_fine = C_state'*f;

figure(1)
clf
if ( length( param1 ) > 1 & length( param2 ) > 1 )
make_logL_contour_plot( log_L, param1, param2, p1_name, p2_name );
else
  semilogx( param1, log_L );
end
title( titlestring );

figure(2)
clf
plot_titration_data( input_data, resnum, conc, pred_fit, sigma_at_each_residue, lane_normalization, conc_fine, pred_fit_fine );
title( titlestring );


for i = 1:size( input_data, 1 )
  input_data_rescale(i,:) = (input_data(i,:) - C_state(1,i))/(C_state(2,i)-C_state(1,i));
  pred_fit_fine_rescale(i,:) = (pred_fit_fine(i,:) - C_state(1,i))/(C_state(2,i)-C_state(1,i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [p_low,p_high] = get_1D_error( param, log_L );

if length( param ) < 2
    p_low  = param;
    p_high = param;
    return;
end

%clf; plot( param, log_L ); pause;

[ log_L_max, min_idx ] = max( log_L )
param( min_idx )
log_L_cutoff = log_L_max - 2.0;

idx = min_idx-1;
while (idx > 0 & log_L(idx) > log_L_cutoff )
    idx = idx - 1;  
end
p_low = param( idx+1 );


idx = min_idx+1;
while (idx <= length(param) & log_L(idx) > log_L_cutoff )
    idx = idx + 1;  
end
p_high = param( idx-1 );







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
