from matplotlib.pylab import *
from rdatkit.datahandlers import RDATFile
from rdatkit.view import VARNA
from rdatkit.secondary_structure import fold
from rdatkit.mapping import MappingData, normalize
from analysis import eigen_reactivities
import sys

rdat = RDATFile()
rdat.load(open(sys.argv[1]))
vals = array(rdat.values.values()[0])
for i in xrange(shape(vals)[0]):
    vals[i,:] = normalize(vals[i,:])
eigenrs = eigen_reactivities(vals)

matshow(vals)
#mshow(vals, cmap=get_cmap('Greys'), vmin=0, vmax=vals.mean(), aspect='auto', interpolation='nearest')
matshow(eigenrs)
#imshow(eigenrs, cmap=get_cmap('Greys'), vmin=eigenrs.min(), vmax=eigenrs.mean(), aspect='auto', interpolation='nearest')
show()
construct = rdat.constructs.values()[0]
for i, e in enumerate(eigenrs[:35]):
    sequence = construct.sequence
    md = MappingData(data=e, seqpos=[s - construct.offset - 1 for s in construct.seqpos])
    print fold(sequence, mapping_data=md)
    structure = fold(sequence, mapping_data=md)[0].dbn
    VARNA.cmd(sequence, structure, 'test_results/eigen_struct%s.png' % i)





