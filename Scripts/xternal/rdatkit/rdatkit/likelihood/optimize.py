import pickle
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.pylab import *
import random
import pdb
from numpy import *
import numpy as np
"""
from cvxmod import *
from cvxmod.atoms import *
from cvxmod.sets import *
"""
from rdatkit import secondary_structure
from rdatkit import settings
from rdatkit import mapping
import scipy.stats as stats
import scipy
import scipy.signal
import scipy.optimize
"""
Next best structure addition (i.e. keep adding structures from the ensemble to optimize)
"""
def by_nbs_addition(sequence, data, database='default.db.dists', use_data=False, nstructs=20):
    db = pickle.load(open('%s/models/%s' % (settings.MAPPING_DATABASE_PATH, database)))
    if use_data:
	structures = secondary_structure.fold(sequence, mapping_data=data, nstructs=nstructs)
    else:
	structures = secondary_structure.fold(sequence, nstructs=nstructs)
    if nstructs > len(structures):
        print 'WARNING: Found %d non trivial structures in the ensemble, not intended %d' % (len(structures), nstructs) 
    struct_probs = []
    struct_lh = []
    for struct in structures:
	probs, lh = struct.likelihood(data, db_obj=db)  
	struct_probs.append(probs)
	struct_lh.append(-log(lh))
    currsol = 0.
    prevsol = float('-Inf')
    numstructs = 1
    while  numstructs <= len(structures):
	"""
	 Solve
	    minimize -sum( ci*log(p(Mi|D)))  over models Mi, given data D, for c
	    0 <= c <= 1 for all i

	    p(Mi|D) are given in struct_lh
	"""
	dim = numstructs
	p = cvxopt.matrix(struct_lh[:dim])
	c = variable(dim, 'c')
	C1 = (c <= 1)
	C2 = (c >= 0)
	C3 = (sum(c) == 1)
	lp = op(min(p.trans()*c), [C1,C2,C3])
	lp.solve()
	print 'Status is %s' % lp.status
	"""
	G = cvxopt.spmatrix(1., range(dim), range(dim))
	h = cvxopt.matrix([1.]*dim)
	A = cvxopt.spmatrix(-1., range(dim), range(dim))
	b = cvxopt.matrix([0.]*dim)
	sol = solvers.lp(c, G, h, A, b)
	opt = array(sol['x'])
	"""
	opt = array(c.value)
	prevsol = currsol
        print 'Solution %s' % opt
        print 'Struct lh %s' % struct_lh[:numstructs]
	currsol = sum(array([opt[i]*struct_lh[i] for i in range(numstructs)]))
	numstructs += 1
    numstructs -= 1
    print 'Finished...'
    print 'Used %d of %d structures' % (numstructs, len(structures))
    return structures[:numstructs], opt, currsol 
    

def by_objective_minimization(sequence, data, t, objfun='norm_l1', structures=[], use_data=False, nstructs=1000):
    print 'Generating structures, this can take a long time...'
    if use_data:
	structures = secondary_structure.sample(sequence, mapping_data=data, nstructs=nstructs, unique=True)
    else:
	if not structures:
	    structures = secondary_structure.sample(sequence, nstructs=nstructs, unique=True)
	else:
	    structures = structures
    if nstructs > len(structures):
        print 'WARNING: Found %d non trivial structures in the ensemble, not intended %d' % (len(structures), nstructs) 
    print 'Preparing inputs for the "Dantzig Selector"'
    A = zeros([len(sequence), len(structures)])
    Aa = np.zeros([len(sequence), len(structures)])
    # CVXMOD having trouble converting to numpy arrays?
    for j, struct in enumerate(structures):
	for i in range(len(sequence)):
	    if struct.dbn[i] == '.':
		A[i,j] = 1.
		Aa[i,j] = 1.
    """
     Solve
	minimize obj_fun( beta )  subject to
	||A'(data-A*beta)||_inf <= (1+(1/t))*sqrt(2*log(len(sequence)))*sigma

	Where sigma is the standard deviation of the data (viewed as a gamma distribution)

    """
    p = len(sequence)
    galpha, gloc, gbeta = stats.gamma.fit(data)
    sigma = sqrt(galpha*(1/gbeta))
    if t != 0:
	l = (1+(1/t))*sqrt(2*log(300))
    else:
	l = 0.
    data[isnan(data)] = 1.
    beta = optvar('b', len(structures))
    Am = param('A', value=A)
    Amt = param('At', value=transpose(A))
    y = param('y', value=matrix(data.tolist()))
    lam = param('lam', value=matrix([l*sigma]))
    con = optvar('con', 1)
    lp = problem()
    print 'Executing the LP'
    if objfun == 'norm_l1':
	lp.constr = [norminf(Amt*(y - Am*beta)) <= lam, beta >= 0 ]
	lp.objective = minimize(norm1(beta))
    if objfun == 'entropy':
	lp.constr = [norminf(Amt*(y - Am*beta)) <= lam, beta >= 0, sum(beta) == 1]
	lp.objective = maximize(entropy(beta))
    if objfun == 'lasso':
	lp.constr = [beta >= 0]
	lp.objective = minimize(norm2(y - Am*beta) + t*norm1(beta))
    print lp
    lp.solve()
    opt = array(beta.value)
    print 'Solution (coefficients) %s' % opt
    error = float(np.linalg.norm(array(data) - np.dot(Aa, opt)))
    predicted = np.dot(Aa, opt)
    print '%d structures have non-zero coefficients' % (np.sum([x > 0 for x in opt[:,0]]))
    return opt, predicted, structures, error

def get_structdistparams(structures, db):
    params = [[(0,0,0)]*len(structures) for i in range(len(structures[0]))]
    structdicts = [s.explode() for s in structures]
    for i in range(len(structures)):
	for k in structdicts[i]:
	    for interval in structdicts[i][k]:
		for idx in interval:
		    params[idx][i] = db[k]
    return params

def conveval(seq1, seq2, seqrange, point):
    conv = scipy.signal.fftconvolve(array(seq1), array(seq2), 'same')
    for i in range(1, len(conv)):
	if point >= seqrange[i-1] and point <= seqrange[i]:
	    return ((conv[i]-conv[i-1])/(seqrange[i]-seqrange[i-1]))*(point - seqrange[i-1]) + conv[i-1]
    return 0.


def structdistconv(distparams, scaleparams, datapoint, convrange, skip=-1):
    dists = []
    for i in range(len(distparams)):
	dists.append(np.vectorize(scipy.stats.gamma(distparams[i][0], loc=distparams[i][1], scale=distparams[i][2]*scaleparams[i]).pdf))
    resrange = dists[0](convrange)
    resrange /= resrange.sum()
    for i in range(1, len(dists)):
	if i != skip:
	    resrange /= resrange.sum()
	    resrange = scipy.signal.fftconvolve(resrange, dists[i](convrange), 'same')
    if datapoint == None:
	return resrange
    for i in range(1, len(convrange)):
	if datapoint >= convrange[i-1] and datapoint <= convrange[i]:
            # return the linear interpolation
	    return ((resrange[i]-resrange[i-1])/(convrange[i]-convrange[i-1]))*(datapoint - convrange[i-1]) + resrange[i-1]
    return 0.

def structdistconv_usingCT(distparams, scaleparams, datapoint, skip=-1):
    means = [0]*len(distparams)
    vars = [0]*len(distparams)
    for i in range(len(distparams)):
        means[i], vars[i] = scipy.stats.gamma.stats(distparams[i][0], loc=distparams[i][1], scale=distparams[i][2]*scaleparams[i], moments='mv')
    res = scipy.stats.norm(loc=array(means).sum(), scale=array(vars).sum()).pdf(datapoint)
    return res
    

def expected_convolution_value(distparams, scaleparams, convrange):
    res = [0]*len(distparams)
    reses = []
    for i in range(len(distparams[1])):
	reses.append([0]*len(distparams))
    print 'enter expected'
    for j in range(len(distparams)):
	dists = []
	for i in range(len(distparams[j])):
	    dists.append(np.vectorize(scipy.stats.gamma(distparams[j][i][0], loc=distparams[j][i][1], scale=distparams[j][i][2]*scaleparams[i]).pdf))
	    reses[i][j] = dot(dists[i](convrange), convrange)
	resrange = dists[0](convrange)
	for i in range(1, len(dists)):
	    resrange /= resrange.sum()
	    resrange = scipy.signal.fftconvolve(resrange, dists[i](convrange), 'same')
	res[j] = dot(resrange, convrange)
	"""
    for i, r in enumerate(reses):
	clf()
	bar(range(len(r)), r)
	savefig('expected_%i.png' % i)
    clf()
    bar(range(len(res)), res)
    savefig('expected_allconv.png')
    pdb.set_trace()
    """
    return res



 
def by_decomvolution(sequence, data, structures, distdb, initcoeffs, learn, tol, maxiterations=100):
    # structdistparams contains the parameters of the distributions of each position in each structure,
    # depending if it's a bulge, helix, etc. structdist[i][j] has the parameters of the distribution of structure j at position i
    structdistparams = get_structdistparams(structures, distdb)
    convrange = np.arange(0., 14, 0.1)
    coeffs = initcoeffs
    iteration = 0
    while True:
	order =  range(len(structures))
	random.shuffle(order)
	currcoeffs = np.copy(coeffs)
	for i in order:
	    print 'In structure %s' % i
	    currcoeffs[i] = coeffs[i]
	    for j in range(len(sequence)):
		shape_i, scale_i, loc_i = structdistparams[j][i]
		fprimevec = np.vectorize(lambda x: (scale_i*(x + loc_i)**(shape_i-1)*(scale_i*currcoeffs[i])**(-shape_i-2)*exp(-(x + loc_i)/(shape_i*currcoeffs[i]))*(x + loc_i - shape_i*scale_i*coeffs[i]))/scipy.special.gamma(shape_i))
		structconv = structdistconv(structdistparams[j], currcoeffs, data[j], convrange)
		convnorm = conveval(fprimevec(convrange), structdistconv(structdistparams[j], currcoeffs, None, convrange, skip=i), convrange, data[j])
		if structconv == 0.:
		    step = 1e-20
		else:
		    step = learn*(1./structconv)*(convnorm)
		currcoeffs[i] += step
		if currcoeffs[i] < 0:
		    currcoeffs[i] = 1e-20
	iteration += 1
        error = np.linalg.norm(coeffs - currcoeffs)
	print 'Finished iteration %s' % iteration
	print 'Current coefficients %s' % currcoeffs
	print 'Iteration distance %s' % error
	coeffs = currcoeffs
	if iteration > maxiterations or error < tol:
	    break
    return coeffs
    


def by_decomvolution_using_CT_approx(sequence, data, structures, distdb, initcoeffs, learn, tol, maxiterations=30, plot_likelihood=False, softmax_reparam=False):
    structdistparams = get_structdistparams(structures, distdb)

    def likelihood_fun(coeffs):
        res = 0
	if softmax_reparam:
	    Z = sum([exp(c) for c in coeffs])
	for i in range(len(sequence)):
	    mu_i = 0
	    sigma_i_sq = 0
	    for j, params in enumerate(structdistparams[i]):
		shape, loc, scale = params
		if softmax_reparam:
		    mu_i += (shape*scale+loc)*exp(coeffs[j])/Z
		    sigma_i_sq += shape*scale**2*(exp(coeffs[j])/Z)**2
		else:
		    mu_i += (shape*scale+loc)*coeffs[j]
		    sigma_i_sq += shape*scale**2*coeffs[j]**2
	    res += 0.5*(log(sigma_i_sq) + ((data[i] - mu_i)/sigma_i_sq)**2)
	return res

    # i is the nucleotide position for which to calculate the likelihood using the gradient
    def likelihood_gradk(coeffs, k, i):
	sigma_i_sq = 0
	mu_i = 0
	if softmax_reparam:
	    Z = sum([exp(c) for c in coeffs])
	for j, params in enumerate(structdistparams[i]):
	    shape, loc, scale = params
	    if softmax_reparam:
		mu_i += (shape*scale+loc)*exp(coeffs[j])/Z
		sigma_i_sq += shape*scale**2*(exp(coeffs[j])/Z)**2
	    else:
		mu_i += (shape*scale+loc)*coeffs[j]
		sigma_i_sq += shape*scale**2*coeffs[j]**2
	    if j == k:
		sigma_ik_sq = shape*scale**2
		mu_ik = shape*scale+loc
	if softmax_reparam:
	    return exp(coeffs[k])/Z*(1 - exp(coeffs[k])/Z)*(sigma_ik_sq/sigma_i_sq - 2*sigma_ik_sq*(data[i] - mu_i)**2/((sigma_i_sq)**3))
	else:
	    return coeffs[k]*sigma_ik_sq/sigma_i_sq + (data[i] - mu_i)/(sigma_i_sq**3)*(-mu_ik*sigma_i_sq - 2*coeffs[k]*sigma_ik_sq*(data[i] - mu_i))

    def likelihood_grad(coeffs):
        res = zeros(len(coeffs))
	for i in range(len(sequence)):
	    for k in range(len(coeffs)):
		res[k] += likelihood_gradk(coeffs, k, i)
	return res

    currcoeffs = array(initcoeffs)
    order =  range(len(sequence))
    iteration = 0
    progcoeffs = [array(initcoeffs)]
    while True:
	random.shuffle(order)
	prevcoeffs = np.copy(currcoeffs)
	#learn = 1./(iteration + 100000)**.5
	for i in order:
	    for k in range(len(structures)):
		currcoeffs[k] -= learn*likelihood_gradk(currcoeffs, k, i)
		if currcoeffs[k] < 0:
		    currcoeffs[k] = random.random()
	iteration += 1
        error = np.linalg.norm(prevcoeffs - currcoeffs)
	print 'Finished iteration %s' % iteration
	print 'Current coefficients %s' % currcoeffs
	print 'Iteration distance %s' % error
	progcoeffs.append(np.copy(currcoeffs))
	if iteration > maxiterations or error < tol:
	    break
    print 'Doing optimization with scipy'
    x0 = np.random.rand(len(initcoeffs))
    print 'Starting x0 for scipy fmin_bfgs %s' % x0
    scipyres = scipy.optimize.fmin_l_bfgs_b(likelihood_fun, x0, fprime=likelihood_grad, iprint=0, pgtol=1e-20, factr=1, bounds=[(0, None) for x in initcoeffs])
    scipyres = scipyres[0]
    print 'Scipy results %s' % scipyres
    if plot_likelihood:
	def likelihood_fun2(c1, c2):
	    return likelihood_fun([c1,c2])
	print 'plotting likelihood function'
	veclik = np.vectorize(likelihood_fun2)
	lim = 7.
	r = arange(0,lim, 0.4)
	pr = array([veclik(r[0], r)])
	for j, i in enumerate(r[1:]):
	    print i
	    pr = np.append(pr, [veclik(i, r)], axis=0)
	"""
	r1, r2 = meshgrid(r,r)
	f = plt.figure()
	ax = Axes3D(f)
	ax.plot_surface(r1, r2, pr)
	f.show()
	pdb.set_trace()
	"""
	clf()
	contour(pr, 10000)
	factor = len(r)/lim
	for i in range(len(progcoeffs)-1):
	    x1 = progcoeffs[i][1]*factor
	    x2 = progcoeffs[i+1][1]*factor
	    y1 = progcoeffs[i][0]*factor
	    y2 = progcoeffs[i+1][0]*factor
	    plot([x1,x2], [y1,y2], 'black')
	scatter([progcoeffs[0][1]*factor, progcoeffs[-1][1]*factor, scipyres[1]*factor], [progcoeffs[0][0]*factor, progcoeffs[-1][0]*factor, scipyres[0]*factor], color=['b', 'r', 'g'])
	#scatter([progcoeffs[0][1]*factor, progcoeffs[-1][1]*factor], [progcoeffs[0][0]*factor, progcoeffs[-1][0]*factor], color=['b', 'r'])
	trange = arange(0,len(r), 3)
	xticks(trange, r[trange])
	yticks(trange, r[trange])
	savefig('likelihoodfun.png')
    print 'Finished optimization'
    print 'Number of iterations %s' % iteration
    print 'Error %s' % error
    print 'Resulting coefficients %s' % currcoeffs
    print 'Likelihood value %s' % likelihood_fun(currcoeffs)
    print 'Gradient values %s' % likelihood_grad(currcoeffs) 
    print 'Resulting coefficients scipy %s' % scipyres
    print 'Likelihood value scipy %s' % likelihood_fun(scipyres)
    print 'Gradient values scipy %s' % likelihood_grad(scipyres) 
    if likelihood_fun(currcoeffs) < likelihood_fun(scipyres):
	print 'Returning stochastic gradient descent coefficients'
	coeffres = currcoeffs
    else:
	print 'Returning bfgs coefficients'
	coeffres = scipyres
    if softmax_reparam:
	Z = sum([exp(c) for c in coeffres])
	for i, c in enumerate(coeffres):
	    coeffres[i] = exp(c)/Z
	print 'Converted coefficients to probabilities %s' % coeffres
    return coeffres

    
def by_decomvolution_using_scipy(sequence, data, structures, distdb, initcoeffs):
    structdistparams = get_structdistparams(structures, distdb)
    convrange = np.arange(0., 14, 0.1)
    def likelihood_fun(coeffs):
        print 'enter likelihood call %s' % coeffs
        res = 0
        for i in range(len(sequence)):
            res -= log(structdistconv_usingCT(structdistparams[i], coeffs, data[i], convrange) + 1e-10)
        print 'exit likelihood call %s' % res
        return res 
    return scipy.optimize.fmin_tnc(likelihood_fun, initcoeffs, approx_grad=True, bounds=[(0, 1) for i in initcoeffs])
    def ge_zero(coeffs):
        for c in coeffs:
            if c < 0:
                return False
        return True
    def sum_one(coeffs):
        s = 0
        for c in coeffs:
            s += c
            if s > 1:
                return False
        return s == 1
    return scipy.optimize.fmin_cobyla(likelihood_fun, initcoeffs, [ge_zero, sum_one])

def by_decomvolution_error_minimization(sequence, data, structures, distb, initcoeffs):
    structdistparams = get_structdistparams(structures, distb)
    convrange = np.arange(0., 14, 0.1)
    data = array(data)
    def error_fun(coeffs):
        return norm(expected_convolution_value(structdistparams, coeffs, convrange) - data)
    return scipy.optimize.fmin_cg(error_fun, initcoeffs)

def by_naive_bayesian(sequence, data, structures, priors, distdb):
    structdistparams = get_structdistparams(structures, distdb)
    posteriors = [0]*len(structures)
    all_shape, all_loc, all_scale = distdb['all']
    for i, s in enumerate(structures):
	probs, probprodb = s.likelihood(data, db_obj=distdb)
	posteriors[i] = prod([probs[j] * priors[i]/ scipy.stats.gamma.pdf(data[j], all_shape, loc=all_loc, scale=all_scale) for j in range(len(sequence))])
    postpartition =  array(posteriors).sum()
    return [p/postpartition for p in posteriors]
