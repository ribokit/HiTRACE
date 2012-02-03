function data_center = centralize( data, refcol, PLOTSTUFF );
% data_center = centralize( data, refcol, PLOTSTUFF );
%
% data      = input data matrix
% refcol    = which column sets the 'starting' normalization (default 1) 
% PLOTSTUFF = make plots of histograms of data ratios, and chosen normalization factors.

niter = 5;
stdev_cut  = 2;

if ~exist( 'refcol' ); refcol = 1; end;
if ~exist( 'PLOTSTUFF' ); PLOTSTUFF = 1; end;

data_center(:,refcol) = data(:,refcol);

if PLOTSTUFF; colorcode = jet( size( data, 2 ) ); r = [-5:0.1:5]; clf; end;

for i  = 2:size( data, 2 )
  logshift = log( abs(data(:,i))./abs(data(:,refcol)) )/log(2);

  % remove anything crazy
  logshift = logshift( find( ~isnan(logshift) & ~isinf( logshift) ) );
  
  if PLOTSTUFF
    h = hist( logshift, r );
    plot( r, h+10*i, 'color', colorcode(i,:) ); hold on
  end
  
  for n = 1:niter
    s = std( logshift );
    m = mean( logshift );
    gp = find(  logshift > (m- stdev_cut*s) & logshift < (m+stdev_cut*s) );
    logshift = logshift( gp );
  end  

  if PLOTSTUFF
    plot( m, interp1(r,h+10*i,m), 'o','color', colorcode(i,:) );  
  end

  data_center(:,i) = data(:,i) / (2^m); 
end
hold off




