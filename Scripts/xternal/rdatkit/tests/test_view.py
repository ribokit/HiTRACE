from rdatkit import view
from rdatkit.secondary_structure import SecondaryStructure
from rdatkit.mapping import MappingData

outfile = open('test_view.html', 'w')
st1 = SecondaryStructure(dbn='....((((....))))....')
st2 = SecondaryStructure(dbn='.....((((..)))).....')
md = MappingData(data=[-0.1]*4+[0.0]*4+[1.0]*4+[0.0]*4+[-1.0]*4)
varna = view.VARNA(['AAAAGGGGUUUUCCCCAAAA'], structures=[st1,st2], mapping_data=[md]) 
outfile.write(varna.render(reference_structure=SecondaryStructure(dbn='(..(((((....)))))..)')))
outfile.close()
