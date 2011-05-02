function [d_sub, bd, alpha] = baseline_matcher( d , d_ref, ymin,ymax)

if ( length( d_ref) == 1 )
  d_ref = d(:,d_ref);
end
if ~exist('ymin')
  ymin= 1;
end
if ~exist('ymax')
  ymax= size(d,1);
end

for k = 1: size( d, 2 );
  fprintf(1,'Baseline-matcher %d\n', k );
  [d_subx, bdx, alphax ] = baseline_matcher_one_profile( d(:,k), d_ref, ...
						  ymin, ymax );
  d_sub(:,k) = d_subx;
  bd(:,k) = bdx;
  alpha( k ) = alphax;
end




function [ d_sub, bd, alpha ] = baseline_matcher_one_profile( d, d_ref, y_min, y_max );

if ~exist( 'ymin' ) 
  ymin = 1;
  ymax = length( d );
end

bd= 0 * d;

% first parameter is A --> weight of smoothness, compared to
%  term favoring as high a baseline as possible.
% Second parameter is B, weight of keeping baseline under signal,
% relative to term favoring as high a baseline as possible.
% Third parameter is C, tries to match profiles to reference lane
% (assumed to be 1);
d_sub = d;
[ d_sub(ymin:ymax), bd(ymin:ymax), alpha ] = baseline_xi( d(ymin:ymax), d_ref(ymin:ymax), ...
					2e3,1, 1e-6 );

%Baseline correction
%d_sub  = ( d  - bd ) / alpha;
bd = d - d_sub;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Baseline Model V1.0
%Yuanxin Xi, IDAV, UC Davis
%Input:
%       b Nx1 the spectrum data
%       A 1x1 smoothing factor
%       B 1x1 negative penalty
%       s 1x1 noise standard deviation
%Output:
%       bd Nx1 baseline
% hacked by Rhiju, Dec. 2009.

function [ gamma_sub, bd, alpha ]=baseline_xi( gamma_in, gamma_ref_in, A,B,C)

%gamma = scale_it( gamma_in );
%gamma_ref = scale_it( gamma_ref_in );
MAX_GAMMA = 100;
gamma = min( gamma_in, MAX_GAMMA);
gamma_ref = min( gamma_ref_in, MAX_GAMMA);

L = length(gamma);

% b is the SIGNAL (gamma in the paper). bd is the estimated baseline (b in the paper).
bd = ones(L,1) * median(gamma); % This was assumed to be ZERO in the paper.
bd0 = gamma;

% current error.
nm = norm( bd - bd0);
nm0 = realmax;

%
% Going to solve:  D * bd = m  
% Note that, relative to paper, D and m have been divided through
% by A. Also, there's a mistake in the expression for M (should be
% +1, not -1).
%
% Also, now solve for a normalization factor "alpha" too -- it will be the
% *last* row/entry in D and m.
%

%M0 = 2  * C * ( gamma - gamma_ref ); 

M0 = 2  * C * gamma;
M0( L+1 ) =   sum( gamma .* gamma_ref );

% initial D matrix
e = ones(L,1);
D0 = spdiags( [ 2*e -8*e (12*e + 2*C/A)  -8*e 2*e], -2:2, L, L );
D0(1,1) = 2; 
D0(L,L) = 2;

D0(2,1) = -4; 
D0(1,2) = -4; 
D0(L,L-1) = -4; 
D0(L-1,L) =-4;

D0(2,2) = 10;
D0(L-1,L-1) = 10;

D0 = D0 * A;

for i = 1:L
  D0( i, L+1 ) =  2 * C * gamma_ref(i);
end

% add entries for  normalization factor (L+1-th entry.)
D0( L+1, L+1 ) = sum( gamma_ref .* gamma_ref );
for i = 1:L
  D0( L+1, i ) = gamma_ref(i);
end

% iteration number
I=0;
alpha = 1.0;

while nm>10 & I<30
  %& nm<nm0;
  I=I+1;

  M=M0;D=D0;bd0=bd;nm0=nm;

  for i=1:L
    if bd(i) > gamma(i)
      M(i) = M(i) + 2 * B * gamma(i);
      D(i,i) = D(i,i) + 2 * B;
    end
  end

  bd=D\M;
  %alpha = bd( L+1 );
  bd = bd(1:L);
  
  % arbitrary convergence check
  nm = norm( bd0 - bd );
end

gamma_sub = (gamma_in - bd) / alpha ;

