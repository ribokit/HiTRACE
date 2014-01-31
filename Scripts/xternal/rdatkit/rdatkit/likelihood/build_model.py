import sys
import pdb
import scipy.stats as stats
import scipy.special
import pickle
from rdatkit.settings import *
from matplotlib.pylab import *
import numpy as np
from scipy import stats, optimize


########## patching scipy

stats.distributions.vonmises.a = -np.pi
stats.distributions.vonmises.b = np.pi


def _fitstart(self, x):
    # method of moment estimator as starting values, not verified
    # with literature
    loc = np.min([x.min(),0])
    a = 4/stats.skew(x)**2
    scale = np.std(x) / np.sqrt(a)
    return (a, loc, scale)

def nnlf_fr(self, thetash, x, frmask):
    # new frozen version
    # - sum (log pdf(x, theta),axis=0)
    #   where theta are the parameters (including loc and scale)
    #
    try:
        if frmask != None:
            theta = frmask.copy()
            theta[np.isnan(frmask)] = thetash
        else:
            theta = thetash
        loc = theta[-2]
        scale = theta[-1]
        args = tuple(theta[:-2])
    except IndexError:
        raise ValueError, "Not enough input arguments."
    if not self._argcheck(*args) or scale <= 0:
        return np.inf
    x = np.array((x-loc) / scale)
    cond0 = (x <= self.a) | (x >= self.b)
    if (np.any(cond0)):
        return np.inf
    else:
        N = len(x)
        #raise ValueError
        return self._nnlf(x, *args) + N*np.log(scale)

def fit_fr(self, data, *args, **kwds):
    loc0, scale0 = map(kwds.get, ['loc', 'scale'],[0.0, 1.0])
    Narg = len(args)
    
    if Narg == 0 and hasattr(self, '_fitstart'):
        x0 = self._fitstart(data)
    elif Narg > self.numargs:
        raise ValueError, "Too many input arguments."
    else:
        args += (1.0,)*(self.numargs-Narg)
        # location and scale are at the end
        x0 = args + (loc0, scale0)
    
    if 'frozen' in kwds:
        frmask = np.array(kwds['frozen'])
        if len(frmask) != self.numargs+2:
            raise ValueError, "Incorrect number of frozen arguments."
        else:
            # keep starting values for not frozen parameters
            x0  = np.array(x0)[np.isnan(frmask)]
    else:
        frmask = None
        
    #print x0
    #print frmask
    return optimize.fmin(self.nnlf_fr, x0,
                args=(np.ravel(data), frmask), disp=0)

stats.distributions.rv_continuous.fit_fr = fit_fr
stats.distributions.rv_continuous.nnlf_fr = nnlf_fr

########## end patching scipy

def KLD(p1, p2):
    return log(scipy.special.gamma(p2[0])) - log(scipy.special.gamma(p1[0])) + p2[0]*(log(p2[2]) - log(p1[2])) + (p1[0]
            - p2[0])*scipy.special.digamma(p1[0]) + p1[0]*(p1[2]-p2[2])/p2[2] + abs(p1[1] - p2[1])
def JSD(p1, p2):
    return (KLD(p1, p2) + KLD(p2, p1))*0.5

db = pickle.load(open(sys.argv[1]))
params = {}
ff = open('oto.txt','w')
for k in db:
    if len(db[k]) > 0:
        vals = array(db[k])
        if k == 'heli':
            p = dists[k].fit_fr(vals[bitwise_and(vals >= 0, vals <= 4)], frozen=[1.0,  0.0, np.nan])
            params[k] = [1, 0, p[0]]
        else:
            p = dists[k].fit_fr(vals[bitwise_and(vals >= 0, vals <= 4)], frozen=[np.nan,  0.0, np.nan])
            params[k] = [p[0], 0, p[1]]
        print 'For %s' % k
        print params[k]
        print 'Mean %s | Variance %s' % (params[k][0]*params[k][2], params[k][0]*(params[k][2]**2))
if len(sys.argv) > 2:
    rnastructfile = open(sys.argv[2] + 'dist.txt', 'w')
    rnastructfile.write('# Gamma distribution parameters: first line is for paired reactivities, second line for unpaired. Order: shape, location and scale\n')
    rnastructfile.write(' '.join([str(x) for x in params['helices']]) + '\n')
    rnastructfile.write(' '.join([str(x) for x in params['unpaired']]) + '\n')
print 'Jensen-Shannon divergence for paired and unpaired distributions: %s' % JSD(params['helices'], params['unpaired'])

pickle.dump(params, open(sys.argv[1][:sys.argv[1].rfind('.')]+'.dists','w'))
