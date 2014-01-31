#!/usr/bin/python

import sys
import os
import pickle
from numpy import *
import argparse
from rdatkit.datahandlers import RDATFile
from rdatkit.secondary_structure import SecondaryStructure
from helpers import normalize
from scipy.stats.stats import scoreatpercentile
import pdb

parser = argparse.ArgumentParser()
parser.add_argument('rdatdir', type=str)
parser.add_argument('outfile', type=str)
parser.add_argument('--normalize', dest='normalize', const=True, default=False, action='store_const')
parser.add_argument('--nooutliers', dest='nooutliers', const=True, default=False, action='store_const')

args = parser.parse_args()

fragtypes = ['all', 'helices', 'interiorloops', 'hairpins', 'dangles', 'bulges',\
        '2wayjunctions', '3wayjunctions', '4wayjunctions', '5wayjunctions', 'unpaired', 'edgepairs', 'internalpairs']
db = {}
dberrors = {}
dbidx = {}
for t in fragtypes:
    db[t] = []
    dberrors[t] = []
    dbidx[t] = {}
for filename in os.listdir(args.rdatdir):
    if not os.path.isdir(args.rdatdir+'/'+filename):
        print filename
    rdat = RDATFile()
    rdat.load(open(args.rdatdir+'/'+filename))
    for cname in rdat.constructs:
        construct = rdat.constructs[cname]
        struct = SecondaryStructure(construct.structure)
        frags = struct.explode()
        for data in construct.data:
            if (('mutation' not in data.annotations) or \
                    ('mutation' in data.annotations and \
                    'WT' in data.annotations['mutation'])):
                if 'modifier' in data.annotations:
                    if args.normalize:
                        normvals = normalize(data.values)
                    else:
                        normvals = data.values
                        iqr = scoreatpercentile(normvals, 75) - scoreatpercentile(normvals, 25)
            for fragtype in frags:
                db['all'].extend(normvals)
                if data.errors:
                    db['all'].extend(data.errors)
                dbidx['all'] = dict([((construct.name, construct.seqpos[i]), v) for i, v in enumerate(normvals)])
                fraglist = frags[fragtype]
                for frag in fraglist:
                    vals = []
                    valerrors = []
                    pos = []
                    for idx in frag:
                        try:
                            iddx = construct.seqpos.index(idx + construct.offset + 1)
                            if ('DMS' in data.annotations['modifier'] and construct.sequence[idx].upper() not in ['A', 'C']) or\
                               ('CMCT' in data.annotations['modifier'] and construct.sequence[idx].upper() not in ['G', 'U']) or\
                                 (args.nooutliers and (normvals[iddx] < 0)):
                                 #(args.nooutliers and (normvals[iddx] > 1.5*iqr or normvals[iddx] < 0)):
                                continue
                            if construct.structure[idx] == '.':
                                db['unpaired'].append(normvals[iddx])
                                dberrors['unpaired'].append(data.errors[iddx])
                                dbidx['unpaired'][(construct.name, idx + construct.offset + 1)] = normvals[iddx]
                            if construct.structure[idx] in (')', '('):
                                db['helices'].append(normvals[iddx])
                                dberrors['helices'].append(data.errors[iddx])
                                dbidx['helices'][(construct.name, idx + construct.offset + 1)] = normvals[iddx]
                                if '.' in (construct.structure[idx-1], construct.structure[idx+1]):
                                    db['edgepairs'].append(normvals[iddx])
                                    dberrors['edgepairs'].append(data.errors[iddx])
                                    dbidx['edgepairs'][(construct.name, idx + construct.offset + 1)] = normvals[iddx]
                                else:
                                    db['internalpairs'].append(normvals[iddx])
                                    dberrors['internalpairs'].append(data.errors[iddx])
                            dbidx['internalpairs'][(construct.name, idx + construct.offset + 1)] = normvals[iddx]
                            val = normvals[iddx]
                            error = data.errors[iddx]
                            if not isnan(val):
                                vals.append(val)
                                valerrors.append(error)
                                pos.append(idx + construct.offset + 1)
                        except ValueError:
                            pass
                    if len(vals) > 0 and fragtype != 'helices':
                        db[fragtype].extend(vals)
                        dberrors[fragtype].extend(valerrors)
                        for i, v in enumerate(vals):
                            dbidx[fragtype][(construct.name, pos[i])] = v

for k, v in db.iteritems():
    f = open(args.outfile + k + '.txt', 'w')
    f.write(','.join([str(x) for x in v]))
pickle.dump(db, open(args.outfile,'w'))
pickle.dump(dbidx, open(args.outfile + '.idx', 'w' ))
pickle.dump(dberrors, open(args.outfile + '.errors', 'w' ))

