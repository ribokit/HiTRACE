function d_correct = attenuation_corrector_v2( d )

d_ratios = 0 * d; %holds fraction of total intensity for each position
d_correct = 0 * d; %holds attenuation - corrected probabilities of reaction at each position

for j = 1:length( d_correct);
   d_ratios(j) =  d(j)/sum( d(1:end));
end;

d_correct(end) = 1;
for j = 1:length(d_correct) - 1;
    d_correct(j) = d_ratios(j) / (sum(d_ratios((j):end)) / sum(d_ratios(1:end )));
end;