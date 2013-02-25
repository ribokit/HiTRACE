function sequence2=RNA2DNA(sequence)

% Converts a RNA sequence to a DNA sequence
% Converts lowercase sequences to an uppercase sequence
% CVL 12/9/10

N_BP = length( sequence );

sequence2=sequence;

for k = 1:N_BP
    c = sequence(k);
    switch c
        case 'U'
            c= 'T';
        case 'u'
            c = 'T';
        case 't'
            c = 'T';
        case 'a'
            c = 'A';
        case 'c'
            c = 'C';
        case 'g'
            c = 'G';
    end
    sequence2( k ) = c;
end
return