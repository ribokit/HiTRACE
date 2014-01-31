import settings
import secondary_structure
from numpy import array

class RNA:
    def __init__(self, sequence=''):
        self.sequence = sequence
    
    def __len__(self):
        return len(self.sequence)

    def bootstrap(self, mapping_data, nboot, nsamples=-1, algorithm='rnastructure', bonus2d=False, replacement=False):
        print 'Starting bootstrap...'
	print 'Folding RNA with complete data'
        if bonus2d:
            mappind_data = array(mapping_data)
        if nsamples < 0:
	    nsamples = len(self.sequence)
        full_bps = secondary_structure.fold(self.sequence, algorithm=algorithm, mapping_data=mapping_data, bonus2d=bonus2d)[0].base_pairs()
	bpdict = dict([(bp, 0) for bp in full_bps])
	for i in range(nboot):
	    print 'Doing bootstrap iteration %s' % i
            if bonus2d:
                grid = indices(mapping_data.shape)
                all_indices = zip(grid[0].ravel(), grid[1].ravel())
                sampled_indices = [choice(all_indices) for x in xrange(len(all_indices))]
                md = zeros(mapping_data.shape)
                for j,k in sampled_indices:
                    md[j,k] += mapping_data[j,k]
            else:
                md = mapping_data.sample(nsamples, replacement=replacement)
	    bps = secondary_structure.fold(self.sequence, algorithm=algorithm, mapping_data=md, bonus2d=bonus2d)[0].base_pairs()
	    for bp in bps:
		if bp in bpdict:
		    bpdict[bp] += 1
		else:
		    bpdict[bp] = 1
        for bp in bpdict:
	    bpdict[bp] *= 100./nboot
	return bpdict

    
