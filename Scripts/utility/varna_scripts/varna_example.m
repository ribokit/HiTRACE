%test

saved_file='ADD_v4/072209_ADD_v4.mat';

importfile(saved_file);
%%
image(area)

sequence='GGAAAGGAAAGGGAAAGAAACGCTTCATATAATCCTAATGATATGGTTTGGGAGTTTCTACCAAGAGCCTTAAACTCTTGATTATGAAGTGAAAACAAAACAAAGAAACAACAACAACAAC';
sequence_actual=sequence(1:end-20);
DBN='....................(((((((((...((((((.........))))))........((((((.......))))))..)))))))))..........';

DATA=area(:,6);
DATA=DATA(end:-1:1);
%%
Z=hSHAPE(DATA);

%%
filename='VARNA/test.html';
varna_fig(filename,sequence_actual,DBN,Z)