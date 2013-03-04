function make_squares( offset, mutpos, native_bp, num_xsel );
%
%  make_squares( offset, mutpos, native_bp, num_xsel );
%

if nargin == 0;  help( mfilename ); return; end;

hold on;

w = 1;
h = 1;
for k = 1: length(mutpos )
  if ( mutpos(k) > 0 )
    x = k-0.5;
    y = 1+num_xsel-(mutpos(k) + offset)-0.5;
    rect_handle = rectangle('Position',[x y w h]);    
  
    gp = find( native_bp(1,:)==mutpos(k));
    if ~isempty(gp)
      for n = gp
	x = k-0.5;
	y = 1+num_xsel-(native_bp(2,n) + offset)-0.5;
	rect_handle = rectangle('Position',[x y w h]);    
	set(rect_handle,'edgecolor',[0.5 0.5 0.5]);
      end
    end
  
  end
end
hold off
