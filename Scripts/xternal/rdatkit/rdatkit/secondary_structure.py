import pdb
import pickle
import os
import tempfile
import scipy.stats as stats
from numpy import *
from random import *
from rdatkit import settings
from rdatkit import mapping


class SecondaryStructure:
    def __init__(self, dbn=''):
        self.dbn = dbn

    def __len__(self):
        return len(self.dbn)

    def base_pairs(self):
        stack = []
        bps = []
        for i, s in enumerate(self.dbn):
            if s == '(':
                stack.append(i)
            if s == ')':
                bps.append((i, stack.pop()))
        return bps

    def base_pair_dict(self):
        stack = []
        bps = {}
        for i, s in enumerate(self.dbn):
            if s == '(':
                stack.append(i)
            if s == ')':
                j = stack.pop()
                bps[j] = i
                bps[i] = j
        return bps


    def helices(self):
        stack = []
        helices = []
        currhelix = []
        for i, s in enumerate(self.dbn):
            if s == '(':
                stack.append(i)
            if s == ')':
                prevbase = stack.pop()
                if len(currhelix) > 0:
                    if currhelix[-1][0] - 1 == prevbase:
                        currhelix.append((prevbase, i))
                    else:
                        helices.append(currhelix)
                        currhelix = [(prevbase, i)]
                else:
                    currhelix.append((prevbase,i))
            if s == '.':
                if len(currhelix) > 0:
                    helices.append(currhelix)
                    currhelix = []
            if len(currhelix) > 0:
                helices.append(currhelix)
        return helices


    def explode(self):
        hstack = [] # for helices
        statestack = [] # for single stranded regions
        jlist = [] # for junctions
        wstack = [] # for junction ways
        jstartstack = []
        fragments = {'helices':[], 'interiorloops':[], 'hairpins':[], 'dangles':[], 'bulges':[]}
        helices = self.helices()
        nhelices = []
        """
        For junction closing
        """
        def _close_junction(junctions, pos):
            junction = []
            ways = 0
            startj = len(junctions)
            for i, jun in enumerate(junctions):
                if jun[0] - 1 == pos:
                    startj = i
                    break
            for i in range(startj, len(junctions)):
                ways += 1
                junction += junctions[i]
            newjunctions = junctions[:startj]
            return junction, ways, newjunctions

        """
        """
        for h in helices:
            nhelices.append([])
            for b in h:
                nhelices[-1].append(b[0])
                nhelices[-1].append(b[1])
        fragments['helices'] = nhelices
        prevstate = '-'
        for i, s in enumerate(self.dbn):
            if s == '(':
                hstack.append(i)
                sstate = ''
                if len(statestack) > 0:
                    sstate = statestack.pop()
                if prevstate == '.' and sstate == '-':
                    fragments['dangles'].append(range(i))
                if (prevstate == '.' or prevstate == ')') and sstate != '-':
                    jlist.append(range(ssregion_start, i))
                    jstartstack.append(ssregion_start-1)
                    wstack.append(1)
                prevstate = '('
            if s == '.':
                if prevstate != '.':
                    ssregion_start = i
                    statestack.append(prevstate)
                prevstate = '.'
            if s == ')':
                prevbase = hstack.pop()
                if prevstate == '.':
                    junction, ways, jlist = _close_junction(jlist, prevbase)
                    if ways > 0:
                        for w in range(ways):
                            jstartstack.pop()
                        if ways == 1:
                            fragments['interiorloops'].append(junction + range(ssregion_start, i))
                        else:
                            key = '%dwayjunctions' % (ways + 1)
                            if key not in fragments:
                                fragments[key] = []
                            fragments[key].append(junction + range(ssregion_start, i))
                    else:
                        sstate = statestack.pop()
                        if sstate == '(':
                            fragments['hairpins'].append(range(ssregion_start, i))
                        if sstate == ')':
                            fragments['bulges'].append(range(ssregion_start, i))
                if prevstate == ')':
                    if len(jstartstack) > 0 and jstartstack[-1] == prevbase:
                        jstartstack.pop()
                        bulge = jlist.pop()
                        fragments['bulges'].append(bulge)
                        wstack.pop()
                        junction_closed = True
                    else:
                        junction, ways, jlist = _close_junction(jlist, prevbase)
                        if ways > 0:
                            for w in range(ways):
                                jstartstack.pop()
                                if ways == 1:
                                    fragments['interiorloops'].append(junction + range(ssregion_start, i))
                                else:
                                    key = '%dwayjunctions' % (ways + 1)
                                    if key not in fragments:
                                        fragments[key] = []
                                    fragments[key].append(junction + range(ssregion_start, i))
                prevstate = ')'
        if self.dbn[-1] == '.':
            fragments['dangles'].append(range(ssregion_start, len(self.dbn)))
        return fragments


    def likelihood(self, mapping_data, database='default.db.dists', db_obj=None):
        frags = self.explode()
        if not db_obj:
            db = pickle.load(open('%s/models/%s' % (settings.MAPPING_DATABASE_PATH, database)))
        else:
            db = db_obj
        data = mapping.normalize(mapping_data)
        probs = array([1.]*len(self.dbn))
        for k in frags:
            if len(frags[k]) > 0 and k in db:
                g = stats.gamma(db[k][0], db[k][1], db[k][2])
                for frag in frags[k]:
                    for i in frag:
                        if g.pdf(data[i]) > 0:
                            probs[i] = g.pdf(data[i])
        return probs, prod(probs)

def to_seqfile(sequence, name='placeholder'):
    seqfile = tempfile.NamedTemporaryFile(delete=False)
    seqfile.write(';\n%(name)s\n%(sequence)s1' % {'name':name, 'sequence':sequence})
    return seqfile

def _prepare_ct_and_seq_files(sequence):
    ctfile = tempfile.NamedTemporaryFile(delete=False)
    ctname = ctfile.name
    ctfile.close()
    seqfile = to_seqfile(sequence)
    seqname = seqfile.name
    seqfile.close()
    return seqname, ctname

def _get_dot_structs(ctname, nstructs, unique=False):
    structs = []
    dbns = []
    for i in range(nstructs):
        dbnfile = tempfile.NamedTemporaryFile(delete=False)
        dbnname = dbnfile.name
        dbnfile.close()
        os.popen(settings.RNA_STRUCTURE_CT2DOT + ' %s %d %s ' % \
             (ctname, i+1, dbnname))
        print(settings.RNA_STRUCTURE_CT2DOT + ' %s %d %s ' % \
             (ctname, i+1, dbnname))
        dbn = open(dbnname).readlines()[-1].strip()
        # Append only non trivial structures
        if '(' in dbn:
            if unique:
                if dbn not in dbns:
                    dbns.append(dbn)
                    structs.append(SecondaryStructure(dbn=dbn))
            else:
                structs.append(SecondaryStructure(dbn=dbn))
    return structs

def _to_ct_file(sequence, struct, filename):
    f = open(filename, 'w')
    if type(struct) != list:
	structs = [struct]
    else:
	structs = struct
    for s in structs:
        f.write('%s bla\n' % len(sequence))
	bps = s.base_pair_dict()
	for i in range(len(sequence)):
	    if i in bps:
		pair = bps[i]+1
	    else:
		pair = 0
	    if i == len(sequence) - 1:
		n = 0
	    else:
		n = i+1
	    f.write('%s %s %s %s %s %s\n' % (i+1, sequence[i], i, i+2, pair, n))
    f.close()

def get_boltzmann_weight(sequence, structure, algorithm='rnastructure'):
    if algorithm == 'rnastructure':
        seqname, ctname = _prepare_ct_and_seq_files(sequence)
        _to_ct_file(sequence, structure, ctname)
        energy = get_energies(ctname)[0]
        PARTCMD = settings.RNA_STRUCTURE_PARTITION + ' %s %s ' % (seqname, '/dev/null')
        print PARTCMD
        ensemble_energy = float(os.popen(PARTCMD).read().split('\n')[-1])
        kT = 0.5905 #TODO Need to generalize/correct this
        weight = exp(-energy/kT)/exp(-ensemble_energy/kT)
    return weight

def mea_structure(sequence, algorithm='rnastructure', nstructs=1, gamma=1.0, opts='', returnct=False):
    if algorithm == 'rnastructure':
        seqname, ctname = _prepare_ct_and_seq_files(sequence)
        CMD = settings.RNA_STRUCTURE + '/exe/MaxExpect -g %s --sequence %s %s %s' % (gamma, seqname, ctname, opts)
        print CMD
        os.popen(CMD)
        structs = _get_dot_structs(ctname, nstructs)
    if returnct:
        return structs, ctname
    else:
        return structs

def fold(sequence, algorithm='rnastructure', mapping_data=[], nstructs=1, fold_opts='', bonus2d=False):
    if algorithm == 'rnastructure':
        seqname, ctname = _prepare_ct_and_seq_files(sequence)
        CMD = settings.RNA_STRUCTURE_FOLD + ' %s %s ' % (seqname, ctname)
        if len(mapping_data) > 0:
            if bonus2d:
                tmp = tempfile.NamedTemporaryFile(delete=False)
                savetxt(tmp, mapping_data)
                tmp.close()
                CMD += '-x %s ' % tmp.name
            else:
                tmp = tempfile.namedtemporaryfile(delete=false)
                tmp.write(str(mapping_data))
                tmp.close()
                cmd += '-sh %s ' % tmp.name
                CMD += fold_opts
        print(CMD)
        os.popen(CMD)
        structs = _get_dot_structs(ctname, nstructs)
    return structs

def partition(sequence, algorithm='rnastructure', mapping_data=None, fold_opts='', bonus2d=False):
    if algorithm == 'rnastructure':
        seqname, ctname = _prepare_ct_and_seq_files(sequence)
        CMD = settings.RNA_STRUCTURE_PARTITION + ' %s %s ' % (seqname, ctname)
        if mapping_data:
            tmp = tempfile.NamedTemporaryFile(delete=False)
	    if bonus2d:
		savetxt(tmp, mapping_data)
		tmp.close()
		CMD += '-x %s ' % tmp.name
	    else:
		tmp.write(str(mapping_data))
		tmp.close()
		CMD += '-sh %s ' % tmp.name
            CMD += fold_opts
        print(CMD)
        os.popen(CMD)
    bppm = loadtxt('bpp.txt')
    for i in range(bppm.shape[0]):
        for j in range(i,bppm.shape[1]):
            if bppm[i,j] != 0:
                bppm[j,i] = bppm[i,j]
    return bppm

def get_structure_energies(sequence, structures):
    energies = []
    ctfile = tempfile.NamedTemporaryFile(delete=False)
    ctname = ctfile.name
    _to_ct_file(sequence, structures, ctname)
    return get_energies(ctname)

def get_energies(ctname):
    EFN2CMD = settings.RNA_STRUCTURE_ENERGY + ' %s /tmp/energy' % ctname
    os.popen(EFN2CMD)
    energyfile = open('/tmp/energy')
    energies = [float(line.split(' ')[-1]) for line in energyfile.readlines()]
    return energies

def sample(sequence, algorithm='rnastructure', mapping_data=None, nstructs=1000, unique=False, energies=False):
    if algorithm == 'rnastructure':
        seqname, ctname = _prepare_ct_and_seq_files(sequence)
        CMD = settings.RNA_STRUCTURE_STOCHASTIC + ' -e %s --sequence %s %s ' % (nstructs, seqname, ctname)
        os.popen(CMD)
        structs = _get_dot_structs(ctname, nstructs, unique=unique)
        energies = get_energies(ctname)
    if energies:
        return structs, energies
    else:
        return structs



def subopt(sequence, algorithm='rnastructure', mapping_data=None, fraction=0.05, nstructs=1000, energies=False):
    if algorithm == 'rnastructure':
        seqname, ctname = _prepare_ct_and_seq_files(sequence)
        CMD = settings.RNA_STRUCTURE_ALLSUB + ' -p %d %s %s ' % (fraction*100, seqname, ctname)
        os.popen(CMD)
        structs = _get_dot_structs(ctname, nstructs)
        energies = get_energies(ctname)
    if energies:
        return structs, energies
    return structs

def random(nstructs, length, nbp):
    structs = []
    for i in xrange(nstructs):
        dbnlist = ['.']*length
        for j in xrange(nbp):
            b2 = randint(0,length-1)
            b1 = randint(0,length-1)
            if b1 <= b2:
                dbnlist[b1] = '('
                dbnlist[b2] = ')'
            else:
                dbnlist[b2] = '('
                dbnlist[b1] = ')'
        structs.append(SecondaryStructure(dbn=''.join(dbnlist)))
    return structs
