function bps = convert_structure_to_bps( structure );

LEFT_BRACKET = [];
bps = [];
for i = 1:length(structure )
  switch structure(i)
   case '('
    LEFT_BRACKET = [LEFT_BRACKET, i];
   case ')'
    bps = [bps; LEFT_BRACKET(end), i];
    LEFT_BRACKET = LEFT_BRACKET( 1 : end-1 );
  end
end