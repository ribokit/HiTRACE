function c2_string = complement( c_string );
%COMPLEMENT(str)
%
%   Returns complement sequence, e.g. AACT > TTGA

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
  c = c_string(k);
  c2 = 'X';
  
  switch c
   case 'A'
    c2 = 'T';
   case 'C'
    c2 = 'G';
   case 'G'
    c2 = 'C';
   case 'U'
    c2= 'A';
   case 'T'
    c2 = 'A';    
  end
  c2_string( k ) = c2;
end
return;
