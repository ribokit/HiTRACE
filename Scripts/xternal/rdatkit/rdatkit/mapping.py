import settings
from scipy.stats import stats
from numpy import *
import random as rand

def normalize(bonuses):
    l = len(bonuses)
    wtdata = array(bonuses)
    if wtdata.min() < 0:
	wtdata -= wtdata.min()
    interquart = stats.scoreatpercentile(wtdata, 75) - stats.scoreatpercentile(wtdata, 25)
    tenperc = stats.scoreatpercentile(wtdata, 90)
    maxcount = 0
    maxav = 0.
    for i in range(l):
	if wtdata[i] >= tenperc:
	    maxav += wtdata[i]
	    maxcount += 1
    maxav /= maxcount
    wtdata = wtdata/maxav
    return wtdata

"""
Returns the maximum likelihood probabilities of adduct formation
and RTase fall-off of the general model described in the paper
of Aviran, et al. "RNA structure characterization from chemical mapping experiments".
X - The data in the plus channel
Y - The data in the minus channel
returns - betas, gammas : the probabilities of adduct formation and RTase fall-off respectively
"""
def maximum_likelihood_probabilities(X, Y):
    betas = [0]*len(X)
    gammas = [0]*len(X)
    sum_X = sum(X)
    sum_Y = sum(Y)
    for i in range(len(X)):
        betas[i] = (X[i]/sum_X - Y[i]/sum_Y)/(1 - sum_Y)
        gammas[i] = Y[i]/sum_Y
    return betas, gammas

def matrix_to_mapping(matrix):
    md = []
    for i in range(shape(matrix)[0]):
            md.append(MappingData(data=matrix[i,:]))
    return md

class MappingData:
    def __init__(self, data=[], seqpos=[], type='', norm=False):
	if seqpos:
	    self._data = [None]*(max(seqpos) + 1)
	    self.seqpos = seqpos
	    for i, pos in enumerate(seqpos):
		self._data[pos] = data[i]
	else:
	    self._data = data
	    self.seqpos = range(len(data))
        if norm:
            self._data = normalize(self._data)
	self.type = type
    
    def data(self):
        return self._data
    
    def load(self, shapefile):
        self.seqpos = []
	mdata = []
        for line in shapefile.readlines():
	    fields = line.strip().split(' ')
	    self.seqpos.append(int(fields[0])-1)
	    mdata.append(float(fields[1]))
	self._data = [None]*(max(self.seqpos) + 1)
	for i, dat in enumerate(mdata):
	    self._data[self.seqpos[i]] = dat
	    
    def __iter__(self):
        return iter(self._data)

    def __len__(self):
        return len(self._data)

    def __getitem__(self, k):
        if k >= len(self._data):
	    return None
        return self._data[k]

    def __str__(self):
        s = ''
        for pos in self.seqpos:
	    if self._data[pos] is not None:
		s += '%d %f\n' % (pos + 1, float(self._data[pos]))
        return s
    
    def sample(self, numsamples, replacement=False):
        if replacement:
	    nseqpos = [0]*numsamples
	    for i in range(numsamples):
		idx = rand.choice(self.seqpos)
		nseqpos[i] = idx
	else:
	    nseqpos = rand.sample(self.seqpos, numsamples)
	ndata = [None]*len(nseqpos)
	for i, pos in enumerate(nseqpos):
	    ndata[i] = self._data[pos]
	return MappingData(data=ndata, seqpos=nseqpos, type=self.type)
	     
