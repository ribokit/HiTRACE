function c2_string = reverse( c_string );
%REVERSE(str)
%
%   Returns reverse sequence, e.g. AACT > TCAA

c_string = upper(c_string);

if iscell( c_string )
  c2_string = {};
  for k = 1: length( c_string )
    c2_string{k} = reverse_complement_string( c_string{k} );
  end
else
  c2_string = reverse_complement_string( c_string );
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
function c2_string = reverse_complement_string( c_string );

for k = 1:length( c_string )
  c = c_string( length(c_string ) - k + 1);
  c2 = c;
  c2_string( k ) = c2;
end
return;
