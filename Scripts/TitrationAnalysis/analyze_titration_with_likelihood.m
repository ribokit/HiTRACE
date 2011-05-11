function [ log_L ] = analyze_titration_with_likelihood( input_data, conc, resnum, K1_conc, nHill, whichres );
%  log_L  = analyze_titration_with_likelihood( input_data, conc, resnum, K1_conc, nHill, whichres );
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
%  nHill   = apparent Hill coefficients to search over. (Default: nHill = [0.0:0.05:4] )
%  whichres = [optional!] just look at this subset of input residues.
%
% (C) Rhiju Das, 2008-2011
%

if ~exist( 'K1_conc' ) | isempty( K1_conc ); K1_conc = 10.^[-3.0 : 0.1 :3.0]; end
if ~exist( 'nHill' )  | isempty( nHill ); nHill = [0.0:0.05:4]; end;
if exist( 'whichres' ) & ~isempty( whichres )
  for k = 1:length(whichres)
    res_to_fit(k) = find( resnum == whichres(k) );
  end
  input_data = input_data( res_to_fit, : );
  resnum = resnum( res_to_fit );
end


if ~exist('colorcode')
  colorcode = [1 0 0; 1 0.5 0; 1 0 1  ; 0.5 0 0    ];
end
if ~exist('markersymbol')
  markersymbol = 'o';
end

best_log_L = -99999999999;
best_K1 = 999;
best_n  = 999;
for i = 1: length( K1_conc)
  for j = 1: length( nHill) % this could be parallelized for speed..
    K1 = K1_conc( i );
    n = nHill( j );
    
    PLOTSTUFF = 0;
    [logL, pred_fit, lane_normalization, sigma_at_each_residue] = ...
	do_new_likelihood_fit( input_data, conc, resnum, ...
			       K1, n );
    
    log_L(i,j) = logL;
    sigma_all(:,i,j) = [sigma_at_each_residue' std( lane_normalization )];
    
    
    fprintf(1,'K1 %8.4f  nHill %8.4f ==>  LogL %8.4f\n',...
	    K1,n, log_L(i,j));
    
    if (log_L(i,j) > best_log_L)
      best_log_L = log_L(i,j);
      best_K1 = K1;
      best_n = n;
    end

%    sigma_at_each_residue'
%    plot_data( input_data, resnum, conc, pred_fit, sigma_at_each_residue );
%    pause;

  end;
end

[logLbest, pred_fit, lane_normalization, sigma_at_each_residue ] =...
    do_new_likelihood_fit( input_data, conc, resnum, ...
			   best_K1, best_n );


% let's try to do a good job of figuring out errors by looking over likelihood
hold on
c = contour( log_L, best_log_L - 2.0 );
K1_low    = interp1( 1:length(K1_conc), K1_conc,  min( c(2,2:end) ) );
K1_high   = interp1( 1:length(K1_conc), K1_conc,  max( c(2,2:end) ) );
n_low    = interp1( 1:length(nHill), nHill,  min( c(1,2:end) ) );
n_high   = interp1( 1:length(nHill), nHill,  max( c(1,2:end) ) );
titlestring = sprintf( 'K1    = %6.4f + %6.4f - %6.4f\n', best_K1, K1_high-best_K1, best_K1 - K1_low );
titlestring = [titlestring, sprintf( 'nHill = %4.2f + %4.2f - %4.2f', best_n, n_high-best_n, best_n - n_low ) ];
fprintf( [titlestring,'\n'] )

figure(1)
clf
make_logL_contour_plot( log_L, K1_conc, nHill );
title( titlestring );

figure(2)
clf
plot_titration_data( input_data, resnum, conc, pred_fit, sigma_at_each_residue, lane_normalization );
title( titlestring );



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
