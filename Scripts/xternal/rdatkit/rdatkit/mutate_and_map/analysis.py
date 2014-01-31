from rdatkit.mapping import *
from matplotlib.pylab import *

def eigen_reactivities(Mmm):
    return svd(Mmm)[2]
