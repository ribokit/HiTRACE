import settings
import sys
from secondary_structure import  *

struct = SecondaryStructure(dbn=sys.argv[1])
print struct.explode()
print len(struct.helices())
print struct.helices()
