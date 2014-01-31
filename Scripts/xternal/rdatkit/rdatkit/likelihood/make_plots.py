from matplotlib.pylab import *
from matplotlib.font_manager import *
import pdb
import pickle
import sys
import scipy.stats as stats
from rdatkit.settings import *

def histOutline(dataIn, *args, **kwargs):
    (histIn, binsIn) = histogram(dataIn, *args, **kwargs)

    stepSize = binsIn[1] - binsIn[0]

    bins = zeros(len(binsIn)*2 + 2, dtype=float)
    data = zeros(len(binsIn)*2 + 2, dtype=float)
    for bb in range(len(binsIn)):
        bins[2*bb + 1] = binsIn[bb]
        bins[2*bb + 2] = binsIn[bb] + stepSize
        if bb < len(histIn):
            data[2*bb + 1] = histIn[bb]
            data[2*bb + 2] = histIn[bb]

    bins[0] = bins[1]
    bins[-1] = bins[-2]
    data[0] = 0
    data[-1] = 0

    return (bins, data)

def KL(p, q):
    res = 0
    for i in range(len(p)):
        res += p[i]*log(p[i]/q[i])
    return res

def JD(p, q):
    m = [(p[i] + q[i])/2 for i in range(len(p))]
    return 0.5*(KL(p,m) + KL(q,m))

figure(2)
title('Combined')
db = pickle.load(open(sys.argv[1]))
params = pickle.load(open(sys.argv[1].replace('.db', '.dists')))
optplot = True
if 'dms' in sys.argv[1]:
    optintercepts = [-1.2,-1,-1,-1,-1, -0.8, -0.8,-0.8,-0.8]
    optslopes = [5.5,4,4.5,5,5.5, 2.5,3,3.5,4]
elif 'cmct' in sys.argv[1]:
    optintercepts = [-1, -1, -0.8, -0.8, -0.8, -0.6, -0.6, -0.6]
    optslopes = [2.5, 3, 1.5, 2, 2.5, 1., 1.5]
elif 'shape' in sys.argv[1]:
    optintercepts = [-1.2,-1.2,-1,-1,-1,-1,-0.8,-0.8,-0.8,-0.8,-0.8,-0.8,-0.8-0.6,-0.6, -0.6]
    optslopes = [3,4,2,2.5,3,3.5,1.5,2,2.5,3,5,5.5,6,1,1.5]
else:
    optplot = False
name = sys.argv[1].replace('.db','')
cl = zip(rand(13), rand(13), rand(13))
all_distobjects = {}
for i, k in enumerate(db):
    if len(db[k]) > 0:
        figure(2)
        if k not in ('all', 'unpaired'):
	    hist(db[k], 100, alpha=0.6, color=cl[i], label=k)
	figure(1)
	clf()
	title(k)
	if dists[k] == stats.gamma:
	    dist = dists[k](params[k][0], loc=params[k][1], scale=params[k][2])
	    #dist = dists[k](params[k][0], scale=params[k][2])
	if dists[k] == stats.expon or dists[k] == stats.norm:
	    dist = dists[k](params[k][0], params[k][1])
        if dists[k] == stats.beta:
	    dist = dists[k](params[k][0], params[k][1], loc=params[k][2], scale=params[k][3])
        all_distobjects[k] = dist
	print 'Plotting %s as %s' % (k, dists[k])
        #plot([dist.pdf(x) for x in frange(min(db[k]), max(db[k]), 0.01)])
        n, bins, patches = hist(db[k], 100, alpha=0.6, color=cl[i], label=k)
        #n, bins, patches = hist(db[k], 100, normed=1, alpha=0.6, color=cl[i])
        #plot(bins, dist.pdf(bins)/dist.pdf(bins).max(), c='r')
	savefig('%s_%s.png' % (name, k))
    else:
	print 'Skipping %s, no data found' % k

prop = FontProperties(size=20)

def set_axis_fontsize():
    fs = 14
    ax = gca()
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(fs)
    ax = gca()
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(fs)
clf()
#title('Paired vs unpaired')
intervals = arange(0, 2, 0.001)
dist = dists['helices'](params['helices'][0], loc=params['helices'][1], scale=params['helices'][2])
bins, npaired = histOutline(db['helices'], 40, range=(0, 2), normed=1)
plot(bins, npaired, color='b', linewidth=2, label='Paired')
#plot(intervals, dist.pdf(intervals), '--k', linewidth=2)
dist = dists['unpaired'](params['unpaired'][0], loc=params['unpaired'][1], scale=params['unpaired'][2])
bins, nunpaired = histOutline(db['unpaired'], 40, range=(0, 2), normed=1)
plot(bins, nunpaired, color='r',linewidth=2,  label='Unpaired')
#plot(intervals, dist.pdf(intervals), '--k', linewidth=3)
#npaired += 1.
#nunpaired += 1
#print('Jensen-Shannon divergence: %s' % JD(npaired/float(sum(npaired)), nunpaired/float(sum(nunpaired))))
ycap = min((99, max(max(npaired), max(nunpaired))))
ylim(0, 7)
xlim(0,2.5)
legend(prop=prop)
set_axis_fontsize()
savefig('%s_paired_vs_unpaired.png' % name)
clf()
title('Edge pairs vs internal pairs')
npaired, bin, patches = hist(db['internalpairs'], 100, alpha=0.6, color='b', label='Internal Pairs')
nunpaired, bin, patches = hist(db['edgepairs'], 100, alpha=0.6, color='r', label='Edge Pairs')
ycap = min((100, max(npaired), max(nunpaired)))
ylim(0, ycap)
xlim(0,2)
legend()
savefig('%s_edge_vs_internal.png' % name)

figure(2)
ylim(0,ycap)
xlim(0,2.3)
legend()
savefig('%s_combined.png' % name)
paireddist = lambda x: dists['helices'](params['helices'][0], loc=params['helices'][1], scale=params['helices'][2]).pdf(x)
print 'paired ' + str(params['helices'])
print 'unpaired ' + str(params['unpaired'])
unpaireddist = lambda x: dists['unpaired'](params['unpaired'][0], loc=params['unpaired'][1], scale=params['unpaired'][2]).pdf(x)
figure(1)
clf()
#title('Probabilstic vs Standard potential')
plot(arange(0,3.5,0.01), [0.5*0.5904*log(unpaireddist(x)/paireddist(x)) for x in arange(0,3.5,0.01)], linewidth=3, color='b',
        label='Probabilistic')
#plot(arange(0,3.5,0.01), [log(x + 1)*2.6 - 0.8 for x in arange(0,3.5,0.01)], color='g', linewidth=2, label='Standard')
if optplot:
    lowbound =  []
    highbound = []
    for i, optparams in enumerate(zip(optintercepts, optslopes)):
        b, m = optparams
        if i == 0:
            lowbound = [log(x + 1)*m + b for x in arange(0,3.5,0.01)]
            highbound = [log(x + 1)*m + b for x in arange(0,3.5,0.01)]
        else:
            curr = [log(x + 1)*m + b for x in arange(0,3.5,0.01)]
            for j, x in enumerate(curr):
                if x > highbound[j]:
                    highbound[j] = x
                if x < lowbound[j]:
                    lowbound[j] = x
    fill_between(arange(0,3.5,0.01), lowbound, highbound, color='g', alpha=0.35)
    plot(arange(0,3.5,0.01), lowbound, color='g', linewidth=2)
    plot(arange(0,3.5,0.01), highbound, color='g', linewidth=2, label='Standard (optimized)')
legend(prop=prop)
ylim(-5, 14)
xlim(0,2)
set_axis_fontsize()
savefig('%s_unpaired_vs_paired_ratio.png' % name)
figure(10)
clf()
#title('Unpaired vs paired distributions')
plot(arange(0,3,0.01), [unpaireddist(x) for x in arange(0,3,0.01)], color='r', linewidth=2, label='Unpaired distribution')
plot(arange(0,3,0.01), [paireddist(x) for x in arange(0,3,0.01)], color='b', linewidth=2, label='Paired distribution')
ylim(0, 5)
legend(prop=prop)
set_axis_fontsize()
savefig('%s_unpaired_vs_paired_dists.png' % name)
