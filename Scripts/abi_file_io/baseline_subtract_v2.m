function [d_sub, bd] = baseline_subtract_v2( d, ymin,ymax, A, B, PLOT_STUFF)
% BASELINE_SUBTRACT_V2:  subtraction of a smooth (but not necessarily constant) baseline.
%
% [d_sub, bd] = baseline_subtract_v2( d, ymin,ymax, A, B, PLOT_STUFF)
%
%  Require input:
%    d = trace (or matrix of traces)
%
%  Optional:
%    ymin, ymax = anchorpoints for subtraction
%    A = parameter weighting how smooth baseline should be [default 2e9]
%    B = parameter weghting how to penalize basleine overshooting input trace. [default 2e4]
%    PLOT_STUFF = setting used by GUI interface to turn off plots.
%
%
% 2010-2011. Optimized & parallelized for nucleic acid capillary electrophoretic traces, R. Das.
%
% Modification of Baseline Model V1.0 Yuanxin Xi, IDAV, UC Davis, originally
%  used for baseline correction of chomatographic traces
% 

d_sub = []; bd = [];
if nargin == 0;  help( mfilename ); return; end;

%Declare boundaries
if ~exist( 'ymin') | ymin == 0;  ymin = 1;end
if ~exist( 'ymax') | ymax == 0;  ymax = size(d,1); end
if ~exist('A') | A == 0;  A = 2e9; end
if ~exist('B') | B == 0;  B = 2e4; end
if ~exist('PLOT_STUFF')  PLOT_STUFF = 1; end


if parallelization_exists()
  %parfor k = 1:size(d,2);
  for k = 1:size(d,2);
      fprintf(1,'Baseline subtracting...%d\n',k);
      [d_sub(:,k),bdx(:,k)] = baseline_subtract_v2_one_profile( d(:,k), ymin,ymax,A,B);
  end 
else
  for k = 1:size(d,2)
      fprintf(1,'Baseline subtracting...%d\n',k);
      [d_sub(:,k),bd(:,k)] = baseline_subtract_v2_one_profile( d(:,k), ymin,ymax,A,B);
  end 
end

if PLOT_STUFF
  subplot(1,2,1);
  scalefactor = 40/mean(mean(d(ymin:ymax,:)));
  image( d(ymin:ymax,:) * scalefactor );
  subplot(1,2,2);
  image( d_sub(ymin:ymax,:) * scalefactor );
  colormap( 1 - gray(100) );
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [d_sub, bd ] = baseline_subtract_v2_one_profile( d, ymin,ymax, ...
						  A,B);

bd=zeros(length(d),1);
%bd(ymin:ymax)=baseline_xi(d(ymin:ymax),2e35,2e29 );

% first parameter is A --> weight of smoothness, compared to
%  term favoring as high a baseline as possible.
% Second parameter is weight of keeping baseline under signal,
% relative to term favoring as high a baseline as possible.
%



bd(ymin:ymax)=baseline_xi(d(ymin:ymax), A,B );

%Baseline correction
d_sub=d - bd;

%Baseline Model V1.0
%Yuanxin Xi, IDAV, UC Davis
%Input:
%       b Nx1 the spectrum data
%       A 1x1 smoothing factor
%       B 1x1 negative penalty
%       s 1x1 noise standard deviation
%Output:
%       bd Nx1 baseline
% 
% new -- filter out negativeoutlier before determining baseline
% 
function bd=baseline_xi(b,A,B,s)

b = filter_negative_outliers( b );

if ( ~exist('s') );   s = 1.0; end

L = length(b);
%Bs = -B/s;  As = -A/s;

% b is the SIGNAL (gamma in the paper). bd is the estimated baseline (b in the paper).
bd = ones(L,1)*median(b); % This was assumed to be ZERO in the paper.
bd0 = b;

% current error.
nm = norm( bd - bd0);
nm0 = realmax;

%
% Going to solve:  D * bd = m  
% Note that, relative to paper, D and m have been divided through
% by A. Also, there's a mistake in the expression for M (should be
% +1, not -1).
%
%M0 = - ones(L,1) / As; 
M0 = s * ones( L, 1 ) / A; 

% initial D matrix
e = ones(L,1);
D0 = spdiags( [ 2*e -8*e 12*e -8*e 2*e], -2:2, L, L );
D0(1,1) = 2; 
D0(L,L) = 2;

D0(2,1) = -4; 
D0(1,2) = -4; 
D0(L,L-1) = -4; 
D0(L-1,L) =-4;

D0(2,2) = 10;
D0(L-1,L-1) = 10;

% iteration number
I=0;

while (nm > 10 | I < 5) & I<30
  %& nm<nm0;
  I=I+1;

  M=M0;D=D0;bd0=bd;nm0=nm;

  for i=1:L
    if bd(i)>b(i)
      %M(i) = M(i) + 2 * Bs * b(i)/As;
      %D(i,i) = D(i,i) + 2 * Bs/As;
      M(i) = M(i) + 2 * (B / A) * b(i);
      D(i,i) = D(i,i) + 2 * ( B / A);
    end
  end

  bd=D\M;
  
  % arbitrary convergence check
  nm = norm( bd0 - bd );
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = filter_negative_outliers( d )
d_sort = sort( smooth( d ) );

N = size( d,1 );
d1 = d_sort( round( N * 0.1 ) );
d3 = d_sort( round( N * 0.9 ) );

%outlier_cutoff = d1 - 1 * abs(d3 - d1 );
outlier_cutoff = d1;% - 0.1*abs(d3-d1);

%plot( d_sort )
%hold on
%plot( [1 N],d1 *[1 1],'r' );
%plot( [1 N],d3 *[1 1],'m' );
%plot( [1 N], outlier_cutoff *[1 1],'k' );
%hold off
%pause;

outlier_points = find( d < outlier_cutoff );
d( outlier_points ) = outlier_cutoff;
