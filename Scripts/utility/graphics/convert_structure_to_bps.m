function bps = convert_structure_to_bps( structure );
%  bps = convert_structure_to_bps( structure );

LEFT_BRACKET = [];

% for pseudoknot (pk) annotation
% pk pairs are in the end of bps
LEFT_BRACKET_PK = [];
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

% for pseudoknot (pk) annotation
% pk pairs are in the end of bps
for i = 1:length(structure )
  switch structure(i)
   case '['
    LEFT_BRACKET_PK = [LEFT_BRACKET_PK, i];
   case ']'
    bps = [bps; LEFT_BRACKET_PK(end), i];
    LEFT_BRACKET_PK = LEFT_BRACKET_PK( 1 : end-1 );
  end
end
