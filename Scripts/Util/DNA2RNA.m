function sequence2 = DNA2RNA( sequence )

% Converts a DNA sequence to an RNA sequence
% Converts lowercase sequences to an uppercase sequence
% CVL 12/9/10

N_BP = length( sequence );

sequence2=sequence;

for k = 1:N_BP
    c = sequence(k);
    switch c
        case 'T'
            c= 'U';
        case 't'
            c = 'U';
        case 'u'
            c = 'U';
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