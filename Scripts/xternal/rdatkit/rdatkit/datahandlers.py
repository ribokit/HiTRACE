"""
datahandlers is a module that contains classes for representing structure mapping file formats in python.
@author Pablo Sanchez Cordero
"""
import ontology
import xlwt
import xlrd
#from Bio import Entrez, Medline
import os
import pdb
from collections import defaultdict
#uncomment if you have scipy installed
#from scipy.cluster.vq import vq, kmeans
"""
Please put your email for Entrez
Entrez.email = ''
"""
def split(s, delims=None):
    for d in delims:
        if d in s:
            return [x for x in s.split(d) if x]

"""
Mapping between RDAT and ISATAB file formats (see sections below)
"""

rdat2isatab_dict = {}


"""
For the RDAT file format for sharing chemical footprinting experiments
"""

class RDATSection:
    pass

class RDATFile:
    def __init__(self):
        self.constructs = defaultdict(list)
        self.traces = defaultdict(list)
        self.values = defaultdict(list)
        self.xsels = defaultdict(list)
        self.data_types = defaultdict(list)
        self.mutpos = defaultdict(list)
	self.errors = defaultdict(list)
	self.annotations = defaultdict(list)
        self.loaded = False

    def append_a_new_data_section( self, current_construct ):
        d = RDATSection()
        d.seqpos = []
        d.errors = []
        d.values = []
        d.xsel = []
        d.trace = []
        d.annotations = {}
        self.constructs[current_construct].data.append(d)


    def load(self, file):
        self.filename = file.name
        current_section = 'general'
        fill_data_types = False
        INIT = False
        self.comments = ''
        while True:
	    line = file.readline()
	    if not line:
		break
	    line = line.strip(' \n')
	    if 'VERSION:' in line:
		self.version = float(line.replace('VERSION:',''))
		continue
	    elif 'VERSION' in line:
		self.version = float(line.replace( 'RDAT_VERSION','').replace('VERSION',''))
		continue
	    if self.version == 0.1:
		if 'COMMENTS:' in line:
		    self.comments = line.replace('COMMENTS:','').strip()
		elif 'ANNOTATION:' in line:
		    if current_section == 'general':
			self.annotations = self.parse_annotations(line.replace('ANNOTATION:', ''))
		    elif current_section == 'construct':
			annotations =  self.parse_annotations(line.replace('ANNOTATION:',''))
			self.constructs[current_construct].annotations =  annotations
			if 'modifier' in annotations:
			    self.data_types[current_construct].append(annotations['modifier'][0])
			    fill_data_types = True
		    elif current_section == 'data':
			annotations = self.parse_annotations(line.replace('ANNOTATION:',''))
			self.constructs[current_construct].data[data_idx].annotations = annotations
			if 'modifier' in annotations:
			    self.data_types[current_construct].append(annotations['modifier'][0])
			if 'mutation' in annotations:
			    try:
				self.mutpos[current_construct][-1] = int(annotations['mutation'][0][1:-1])
			    except ValueError:
				pass
		    else:
			print 'Attribute :'+line+' does not belong to a valid section'

		elif 'CONSTRUCT' in line:
		    current_section = 'construct'
		    if fill_data_types:
                        self.data_types[current_construct] = [self.data_types[current_construct][0]]*len(self.values[current_construct])
		    	fill_data_types = False
		    current_construct = file.readline().strip().replace('NAME:','').strip()
		    data_idx = -1
		    self.constructs[current_construct] = RDATSection()
		    self.constructs[current_construct].name = current_construct
		    self.constructs[current_construct].data = []
		    self.constructs[current_construct].seqpos = []
		    self.constructs[current_construct].xsel = []
		    self.traces[current_construct] = []
		    self.values[current_construct] = []
		    self.xsels[current_construct] = []
		    self.data_types[current_construct] = []
		    self.mutpos[current_construct] = []
		    self.constructs[current_construct].structure = ''
		elif 'SEQUENCE:' in line:
		    self.constructs[current_construct].sequence = line.replace('SEQUENCE:', '').strip()
		elif 'STRUCTURE:' in line:
		    self.constructs[current_construct].structure = line.replace('STRUCTURE:', '').strip()
		elif 'WELLS:' in line:
		    self.constructs[current_construct].wells = line.replace('WELLS:', '').strip().split(',')
		elif 'OFFSET:' in line:
		    self.constructs[current_construct].offset = int(line.replace('OFFSET:',''))
		elif 'DATA' in line:
		    current_section = 'data'
		    data_idx += 1
		    d = RDATSection()
		    d.xsel = []
		    d.seqpos = []
		    self.mutpos[current_construct].append('WT')
		    self.constructs[current_construct].data.append(d)
		elif 'SEQPOS:' in line:
		    if current_section == 'construct':
			    self.constructs[current_construct].seqpos= [int(x) for x in line.replace('SEQPOS:','').strip(' ,').split(',')]
		    else:
			    self.constructs[current_construct].data[data_idx].seqpos = [int(x) for x in line.replace('SEQPOS:','').strip(' ,').split(',')]
		elif 'VALUES:' in line:
		    self.constructs[current_construct].data[data_idx].values = [float(x) for x in line.replace('VALUES:','').strip(' ,').split(',')]
		    self.values[current_construct].append(self.constructs[current_construct].data[data_idx].values)
		elif 'TRACE:' in line:
		    self.constructs[current_construct].data[data_idx].trace = [float(x) for x in line.replace('TRACE:','').strip(' ,').split(',')]
		    self.traces[current_construct].append(self.constructs[current_construct].data[data_idx].trace)
		elif 'XSEL:' in line:
		    if current_section == 'construct':
			    self.constructs[current_construct].xsel = [float(x) for x in line.replace('XSEL:','').strip(' ,').split(',')]
		    else:
			    self.constructs[current_construct].data[data_idx].xsel = [float(x) for x in line.replace('XSEL:','').strip(' ,').split(',')]
			    self.xsels[current_construct].append(self.constructs[current_construct].data[data_idx].xsel)
		else:
		    print 'Invalid section: '+line
	    elif self.version >= 0.2:
		if 'COMMENT' in line:
		    self.comments += line.replace('COMMENTS','').replace('COMMENT ','') + '\n'
		elif 'ANNOTATION' in line and not 'ANNOTATION_DATA' in line:
		    self.annotations = self.parse_annotations(split(line.replace('ANNOTATION', ''), delims='\t '))
		elif 'CONSTRUCT' in line or 'NAME' in line:
		    current_section = 'construct'
                    if 'CONSTRUCT' in line: line = file.readline().strip() # Advance to 'NAME' line.
		    #if fill_data_types:
		    #	self.data_types[current_construct] = [self.data_types[current_construct][0]]*len(self.values[current_construct])
		    #	fill_data_types = False
		    current_construct = line.replace('NAME','').strip()
		    data_idx = -1
		    self.constructs[current_construct] = RDATSection()
		    self.constructs[current_construct].name = current_construct
		    self.constructs[current_construct].data = []
		    self.constructs[current_construct].seqpos = []
		    self.constructs[current_construct].xsel = []
		    self.constructs[current_construct].annotations = {}
		    self.traces[current_construct] = []
		    self.values[current_construct] = []
		    self.xsels[current_construct] = []
		    self.data_types[current_construct] = []
		    self.mutpos[current_construct] = []
		    self.errors[current_construct] = []
		    self.constructs[current_construct].structure = ''
		    self.constructs[current_construct].sequence = ''
		    self.constructs[current_construct].structures =  defaultdict(int)
		    self.constructs[current_construct].sequences =  defaultdict(int)
		elif 'SEQUENCE' in line:
		    if ':' in line:
			seqidx, seq = line.replace('SEQUENCE:','').strip().split()
			self.constructs[current_construct].sequences[int(seqidx)] = seq.strip()
		    else:
			self.constructs[current_construct].sequence = line.replace('SEQUENCE', '').strip()
                elif 'STRUCTURE' in line:
		    if 'STRUCTURE:' in line:
			structidx, struct = split(line.replace('STRUCTURE:','').strip(), delims='\t ')
			self.constructs[current_construct].structures[int(structidx)] = struct.strip()
		    else:
			self.constructs[current_construct].structure = line.replace('STRUCTURE', '').strip()
		elif 'OFFSET' in line:
		    self.constructs[current_construct].offset = int(line.replace('OFFSET',''))
		elif 'DATA_TYPE' in line:
		    current_section = 'data'
		    self.data_types[current_construct] = split(line.replace('DATA_TYPE', '').strip(), delims='\t ')
		elif 'SEQPOS' in line:
                    if self.version >= 0.32:
                        self.constructs[current_construct].seqpos= [int(x[1:]) for x in split(line.replace('SEQPOS','').strip(), delims='\t ')]
                    else:
                        self.constructs[current_construct].seqpos= [int(x) for x in split(line.replace('SEQPOS','').strip(), delims='\t ')]
		elif 'MUTPOS' in line:
		    self.mutpos[current_construct] = [x.strip() for x in split(line.replace('MUTPOS','').strip(), delims='\t ')]
		elif 'ANNOTATION_DATA' in line:
		    if self.version >= 0.23:
			fields = split(line.replace('ANNOTATION_DATA:', '').strip(), delims='\t ')
		    else:
			fields = split(line.replace('ANNOTATION_DATA ', '').strip(), delims='\t ')
		    data_idx = int(fields[0])-1
		    annotations = self.parse_annotations(fields[1:])
                    for l in xrange(data_idx - len(self.constructs[current_construct].data) + 1):
                        self.append_a_new_data_section( current_construct )
		    self.constructs[current_construct].data[data_idx].annotations = annotations
		    if 'mutation' in annotations:
			try:
			    if len(self.mutpos[current_construct]) > 0:
				self.mutpos[current_construct][-1] = int(annotations['mutation'][0][1:-1])
		            else:
				self.mutpos[current_construct].append(int(annotations['mutation'][0][1:-1]))
			except ValueError:
			    pass
		elif 'AREA_PEAK ' in line or 'REACTIVITY:' in line:
		    if 'AREA_PEAK ' in line:
			fields = split(line.replace('AREA_PEAK ', '').strip('\n ,'), delims='\t ')
		    if 'REACTIVITY:' in line: # Means we are in version >= 0.23
			fields = split(line.replace('REACTIVITY:', '').strip('\n ,'), delims='\t ')
		    data_idx = int(fields[0])-1
		    peaks = [float(x) for x in fields[1:]]
                    if ( data_idx >= len(self.constructs[current_construct].data) ):
                        self.append_a_new_data_section( current_construct )
		    self.constructs[current_construct].data[data_idx].values = peaks
		    self.values[current_construct].append(self.constructs[current_construct].data[data_idx].values)
		elif 'AREA_PEAK_ERROR' in line or 'REACTIVITY_ERROR:' in line:
		    if 'AREA_PEAK_ERROR ' in line:
			fields = split(line.replace('AREA_PEAK_ERROR ', '').strip('\n ,'), delims='\t ')
		    if 'REACTIVITY_ERROR:' in line: #Means we are in version >= 0.23
			fields = split(line.replace('REACTIVITY_ERROR:', '').strip('\n ,'), delims='\t ')
		    data_idx = int(fields[0])-1
		    errors = [float(x) for x in fields[1:]]
		    self.constructs[current_construct].data[data_idx].errors = errors
		    self.errors[current_construct].append(self.constructs[current_construct].data[data_idx].errors)
		elif 'TRACE' in line:
		    fields = split(line.replace(':', ' ').replace('TRACE', '').strip('\n ,'),delims='\t ')
		    data_idx = int(fields[0])-1
		    trace = [float(x) for x in fields[1:]]
		    self.constructs[current_construct].data[data_idx].trace = trace
		    self.traces[current_construct].append(self.constructs[current_construct].data[data_idx].trace)
		elif 'XSEL_REFINE' in line:
		    fields = split(line.replace(':', ' ').replace('XSEL_REFINE', '').strip('\n ,'),delims='\t ')
		    data_idx = int(fields[0])-1
		    xsel = [float(x) for x in fields[1:]]
		    self.constructs[current_construct].data[data_idx].xsel = xsel
		    self.xsels[current_construct].append(self.constructs[current_construct].data[data_idx].xsel)
		elif 'XSEL' in line:
		    self.constructs[current_construct].xsel = [float(x) for x in split(line.replace(':', ' ').replace('XSEL','').strip('\n ,'),delims='\t ')]
		else:
                    if line.strip():
			    print 'Invalid section: '+line
        else:
	    print 'Unknown version %s!' % self.version
        if fill_data_types:
            self.data_types[current_construct] = [self.data_types[current_construct][0]]*len(self.values[current_construct])
        self.loaded = True


    def parse_annotations(self, s):
        d = {}
	if self.version == 0.1:
	    token = ';'
            s = s.split(',')
	else:
	    token = ':'
	for item in s:
            if item:
		pair = item.split(token)
		if pair[0].strip() in d:
		    d[pair[0].strip()].append(':'.join(pair[1:]))
		else:
		    d[pair[0].strip()] = [':'.join(pair[1:])]
	return d

    def annotation_str(self, a, delim):
        s = ''
        for k in a:
	    for i in range(len(a[k])):
	        print a[k]
	        s += k + ':' + a[k][i] + delim
	return s

    def save(self, filename, version=None, delim='\t'):
        if not version:
	    if not self.version:
		version = 0.3
	    else:
		version = self.version
        if not self.loaded:
	    print 'No data to save...'
	else:
	    if version == 0.1:
		file = open(filename, 'w')
		file.write('VERSION: '+str(self.version)+'\n')
		file.write('COMMENTS: '+str(self.comments)+'\n')
		file.write('ANNOTATION: '+self.annotation_str(self.annotations, delim)+'\n')
		for name in self.constructs:
		    construct = self.constructs[name]
		    file.write('CONSTRUCT\n')
		    file.write('NAME: '+name+'\n')
		    file.write('SEQUENCE: '+construct.sequence+'\n')
		    file.write('WELLS: '+','.join([x for x in construct.wells])+'\n')
		    file.write('OFFSET: '+str(construct.offset)+'\n')
		    file.write('ANNOTATION: '+self.annotation_str(construct.annotations, delim)+'\n')
		    for d in construct.data:
			file.write('DATA\n')
			file.write('SEQPOS: '+','.join([str(x) for x in d.seqpos])+'\n')
			file.write('ANNOTATION: '+self.annotation_str(d.annotations, delim)+'\n')
			file.write('VALUES: '+','.join([str(x) for x in d.values])+'\n')
            elif version == 0.2 or version == 0.21 :
	        file = open(filename, 'w')
		file.write('VERSION '+str(version)+'\n')
		file.write('COMMENTS '+str(self.comments)+'\n')
		for name in self.constructs:
		    construct = self.constructs[name]
		    file.write('CONSTRUCT\n')
		    file.write('NAME '+name+'\n')
		    file.write('SEQUENCE '+construct.sequence+'\n')
		    file.write('STRUCTURE '+construct.structure+'\n')
		    file.write('OFFSET '+str(construct.offset)+'\n')
		    if construct.annotations:
			file.write('ANNOTATION '+self.annotation_str(construct.annotations, delim)+'\n')
		    file.write('MUTPOS '+delim.join([str(x) for x in self.mutpos[name]])+'\n')
		    file.write('SEQPOS '+delim.join([str(x) for x in construct.seqpos])+'\n')
		    file.write('DATA_TYPE '+' '.join(self.data_types[name])+'\n')
		    for i, d in enumerate(construct.data):
			file.write('ANNOTATION_DATA %s %s\n' % (i+1, self.annotation_str(d.annotations, delim)))
		    for i,row in enumerate(self.values[name]):
			file.write('AREA_PEAK %s %s\n' % (i+1, delim.join([str(x) for x in row])))
		    for i,row in enumerate(self.traces[name]):
			file.write('TRACE %s %s\n' % (i+1, delim.join([str(x) for x in row])))
		    if self.errors:
			for i,row in enumerate(self.errors[name]):
			    file.write('AREA_PEAK_ERRORS %s %s\n' % (i+1, delim.join([str(x) for x in row])))
		    if construct.xsel:
			file.write('XSEL %s\n' % ' '.join([str(x) for x in construct.xsel]))
		    for i,row in enumerate(self.xsels[name]):
			file.write('XSEL_REFINE %s %s\n' % (i+1, delim.join([str(x) for x in row])))
            elif version >= 0.24 :
	        file = open(filename, 'w')
		file.write('VERSION%s' % (delim) + str(version)+'\n')
		file.write('COMMENTS%s' % (delim) + str(self.comments)+'\n')
		for name in self.constructs:
		    construct = self.constructs[name]
		    file.write('NAME%s' % (delim) + name+'\n')
		    if construct.sequence:
			file.write('SEQUENCE%s' % (delim) + construct.sequence+'\n')
		    for k, v in construct.sequences.iteritems():
			file.write('SEQUENCE:%s%s%s\n' % (k, delim, v))
		    if construct.structure:
			file.write('STRUCTURE%s' % (delim) + construct.structure+'\n')
		    for k, v in construct.structures.iteritems():
			file.write('STRUCTURE:%s%s%s\n' % (k, delim, v))
		    file.write('OFFSET%s' % (delim) + str(construct.offset)+'\n')
		    if construct.annotations:
			file.write('ANNOTATION%s' % (delim) + self.annotation_str(construct.annotations, delim)+'\n')
                    if version >= 0.32:
                        if name in self.mutpos:
                            file.write('MUTPOS%s' % (delim) + delim.join([str(x) for x in self.mutpos[name]])+'\n')
                        else:
                            file.write('MUTPOS%s' % (delim) + 'WT '*len(construct.data) + '\n')
                    if version >= 0.32:
                        file.write('SEQPOS%s' % (delim) + delim.join([construct.sequence[x - construct.offset] + str(x+1) for x in construct.seqpos])+'\n')
                    else:
                        file.write('SEQPOS%s' % (delim) + delim.join([str(x+1) for x in construct.seqpos])+'\n')
		    for i, d in enumerate(construct.data):
			file.write('ANNOTATION_DATA:%s%s%s\n' % (i+1, delim, self.annotation_str(d.annotations, delim)))
		    if name in self.values:
			for i,row in enumerate(self.values[name]):
			    file.write('REACTIVITY:%s%s%s\n' % (i+1, delim, delim.join([str(x) for x in row])))
		    if name in self.traces:
			for i,row in enumerate(self.traces[name]):
			    if len(row) > 0:
				file.write('TRACE:%s%s%s\n' % (i+1, delim, delim.join([str(x) for x in row])))
		    if name in self.errors:
			for i,row in enumerate(self.errors[name]):
			    if len(row) > 0:
				file.write('REACTIVITY_ERRORS:%s%s%s\n' % (i+1, delim, delim.join([str(x) for x in row])))
		    if construct.xsel:
			file.write('XSEL%s%s\n' % (delim, delim.join([str(x) for x in construct.xsel])))
		    if name in self.xsels:
			for i,row in enumerate(self.xsels[name]):
			    if len(row) > 0:
				file.write('XSEL_REFINE:%s%s%s\n' % (i+1, delim, delim.join([str(x) for x in row])))

	    else:
		print 'Wrong version number %s' % version
    
    def validate(self):
        messages = []
        for name in self.constructs:
	    c = self.constructs[name]
	    if len(name) == 0:
		messages.append( 'WARNING! Must give a name!')
	    if len(c.sequence) == 0:
		messages.append( 'WARNING! Must supply sequence!')
	    if 'T' in c.sequence:
		messages.append( 'WARNING! Warning: you have a T instead of a U in the sequence!!')
	    if min(c.seqpos) - c.offset < 1:
		messages.append( 'WARNING! Offset/seqpos does not look right -- at least one index is too low for sequence')
	    if max(c.seqpos) - c.offset > len(c.sequence):
		messages.append( 'WARNING! Offset/seqpos does not look right -- at least one index is too high for sequence') 
	    if len(c.data[0].values) != len(c.seqpos):
                messages.append( 'WARNING! Number of bands in area_peak [%s] does not match len of seqpos [%s]' % (len(c.data[0].values), len(c.seqpos))) 
            for i, d in enumerate(c.data):
		if not 'annotations' in d.__dict__:
		    messages.append( 'WARNING! Data for index %s has no annotations' % i)
		if not 'values' in d.__dict__:
		    messages.append( 'WARNING! Data for index %s has no values for area peaks' % i)
		if not 'trace' in d.__dict__:
		    messages.append( 'WARNING! Data for index %s has no trace' % i)
		if len(self.xsels) > 0 and not 'xsel' in d.__dict__:
		    messages.append( 'WARNING! Data for index %s has no xsel refine' % i)
		if 'xsel' in c.__dict__:
		    if len(c.xsel) != len(d.values):
                        messages.append( 'WARNING! Number of bands in construct xsel [%s] does not match number of bands in values area peak [%s] of data indexed %s' % (len(c.xsel), len(d.values), i ))
	        if 'xsel' in d.__dict__: 
		    if len(d.xsel) != 0 and len(d.xsel) != len(d.values):
                        messages.append( 'WARNING! Number of bands in xsel indexed %s  [%s] does not match number of bands in values area peak [%s]' % (i, len(d.xsel), len(d.values) ))
	    return messages
 
		
    def toISATAB(self):
        def parse_concentration(s):
	    for i, ch in enumerate(s):
		if not ch in [str(x) for x in range(10) + ['.']]:
		    idx = i
		    break
            return s[:idx], s[idx:]
	isatabfile = ISATABFile()
	j = 0
	general_factors = []
	protocols = set()
	general_protocol = ''
	chemicals = set()
	isatabfile.investigation_dict['Study File Name'].append(self.filename + ' (check entry ' + self.filename[:self.filename.find('.')] + ' at http://rmdb.stanford.edu for details)')
	"""
        if 'pmid' in self.annotations:
	    pmid = self.annotations['pmid'][0]
	    h = Entrez.efetch(db='pubmed', id=[pmid], rettype='medline', retmode='text')
	    records = Medline.parse(h)
	    for r in records:
		record = r
	    doi = ''
	    for item in record.get('AID', '?'):
		if '[doi]' in item:
		    doi = item.replace('[doi]','')
	    isatabfile.investigation_dict['Study Title'].append(record.get('TI', '?'))
	    isatabfile.investigation_dict['Study Publication DOI'].append(doi)
	    isatabfile.investigation_dict['Study Publication Title'].append(record.get('TI', '?'))
	    isatabfile.investigation_dict['Study Publication Author list'].append(','.join(record.get('AU', '?')))
	    isatabfile.investigation_dict['Study Description'].append(record.get('AB', '?'))
	    isatabfile.investigation_dict['Study Public Release Date'].append(record.get('DP', '?'))
	    isatabfile.investigation_dict['Study PubMed ID'].append(pmid)
	    isatabfile.investigation_dict['Study Publication Status'].append('indexed in pubmed')
        """
	for k in ['chemical', 'salt', 'buffer', 'temperature']:
	    if k in self.annotations:
		isatabfile.investigation_dict['Study Factor Name'].append('%s %s (constant for all assays)' % (k,self.annotations[k]))
		if k != 'temperature':
		    isatabfile.investigation_dict['Study Factor Type'].append('Compound')
		else:
		    isatabfile.investigation_dict['Study Factor Type'].append('Temperature')
	if 'technology' in self.annotations:
	    tech = self.annotations['technology'][0]
	else:
	    tech = 'capillary electrophoresis'
	if 'modifier' in self.annotations:
	    if self.annotations['modifier'][0] in ontology.modifier_protocol:
		general_protocol = ontology.modifier_protocol[self.annotations['modifier'][0]]
            else:
		general_protocol = self.annotations['modifier'][0]
	    protocols.add(general_protocol)

	"""
	for k in ['chemical', 'salt', 'buffer']:
	    if k in self.annotations:
		name, concentration = self.annotations[k][0].split(':')
		isatabfile.investigation_dict['Study Factor Name'].append(name)
		isatabfile.investigation_dict['Study Factor Name'].append(name + ' concentration')
		term = ontology.chemicals[name]
		concentration, units = parse_concentration(concentration)
		general_factors.append(['Factor Value[%s]' % name, name])
		general_factors.append(['Term Source REF[%s]' % name, term.split(':')[0]])
		general_factors.append(['Term Accession Number[%s]' % name, term])
		general_factors.append(['Factor Value[%s concentration]' % name, name])
		general_factors.append(['Unit[%s concentration]' % name, units])
		isatabfile.assays_keys.append('Factor Value[%s]' % name)
		isatabfile.assays_keys.append('Term Source REF[%s]' % name)
		isatabfile.assays_keys.append('Term Accession Number[%s]' % name)
		isatabfile.assays_keys.append('Factor Value[%s concentration]' % name)
		isatabfile.assays_keys.append('Unit[%s concentration]' % name)
	"""
	for cname in self.constructs:
	    construct = self.constructs[cname] 
	    for i, d in enumerate(construct.data):
		name = cname.replace(' ','-')
		seq = ''
		for j in construct.seqpos:
		    seq += construct.sequence[j - construct.offset - 1]
		if 'mutation' in d.annotations:
		    for j in range(len(d.annotations['mutation'])):
			mutlabel = d.annotations['mutation'][j].strip('\t')
			name = name + '_' + mutlabel
			if mutlabel != 'WT':
			    try:
				index = int(mutlabel[1:-1]) - construct.offset
				seq = seq[:index-1]+mutlabel[-1]+seq[index:]
			    except ValueError:
			        # The mutation label is not in standard format, default to normal sequence and make a note
				seq += ' Note, mutation=%s' % mutlabel
		else:
		    name = name + '_WT'
		idname = name + '_' + str(i+1)
		isatabfile.assays_dict['Assay Name'].append(idname)
		isatabfile.sample_id_name_map[idname] = name
		isatabfile.data[idname] = d.values
		isatabfile.data_id_order.append(idname)
		isatabfile.assays_dict['Source Name'].append(name)
		isatabfile.assays_dict['Characteristics[Nucleotide Sequence]'].append(seq)
		isatabfile.assays_dict['Characteristics[Nucleotide Type]'].append('RNA')
		if 'production' in d.annotations:
		    isatabfile.assays_dict['Characteristics[RNA Production]'].append(d.annotations['production'][0].replace('-',' '))
		else:
		    isatabfile.assays_dict['Characteristics[RNA Production]'].append('in vitro synthesis')
		if 'modifier' in d.annotations:
		    modifier = d.annotations['modifier'][0]
		    if modifier in ontology.modifier_protocol:
			isatabfile.assays_dict['Protocol REF'].append(ontology.modifier_protocol[d.annotations['modifier'][0]].replace('-',' '))
			protocol = ontology.modifier_protocol[d.annotations['modifier'][0]]
		    else:
			isatabfile.assays_dict['Protocol REF'].append(modifier)
			protocol = modifier
		    protocols.add(protocol)
		elif 'modifier' in self.annotations:
                    modifier = self.annotations['modifier'][0]
		    if modifier in ontology.modifier_protocol:
			isatabfile.assays_dict['Protocol REF'].append(ontology.modifier_protocol[self.annotations['modifier'][0]].replace('-',' '))
			protocol = ontology.modifier_protocol[self.annotations['modifier'][0]]
	            else:
			isatabfile.assays_dict['Protocol REF'].append(modifier)
			protocol = modifier
		    protocols.add(protocol)
		elif general_protocol:
		    isatabfile.assays_dict['Protocol REF'].append(general_protocol)
		    protocol = general_protocol
		else:
		    isatabfile.assays_dict['Protocol REF'].append('')
		    protocol = ''
		isatabfile.assays_dict['Parameter Value[Data Starts at Sequence Position]'].append(str(min(construct.seqpos)-construct.offset))
		isatabfile.assays_dict['Raw Data File'].append('datamatrix.txt')
		if 'performer' in d.annotations:
		    isatabfile.assays_dict['Performer'].append(d.annotations['performer'][0].replace('-',' '))
		elif 'performer' in self.annotations:
		    isatabfile.assays_dict['Performer'].append(self.annotations['performer'][0].replace('-',' '))
		else:
		    isatabfile.assays_dict['Performer'].append('')
		if 'date' in d.annotations:
		    isatabfile.assays_dict['Date'].append(d.annotations['date'][0])
		elif 'date' in self.annotations: 
		    isatabfile.assays_dict['Date'].append(self.annotations['date'][0])
		else:
		    isatabfile.assays_dict['Date'].append('')
		isatabfile.assays_dict['Term Source REF'].append('OBI')
		if protocol in ontology.protocols:
		    isatabfile.assays_dict['Term Accession Number'].append(ontology.protocols[protocol])
		else:
		    isatabfile.assays_dict['Term Accession Number'].append(protocol)
		for k in ['chemical', 'folding-salt', 'buffer', 'salt']:
		    if k in d.annotations:
			chemical, concentration = d.annotations[k][0].split(':')
			concentration, units = parse_concentration(concentration)
                        if chemical in ontology.chemicals:
		            term = ontology.chemicals[chemical]
			else:
			    term = chemical
			if not k in isatabfile.assays_factors:
			    isatabfile.assays_factors[k] = {}
			    isatabfile.assays_factors[k]['value'] = []
			    isatabfile.assays_factors[k]['ref'] = []
			    isatabfile.assays_factors[k]['concentration'] = []
			    isatabfile.assays_factors[k]['unit'] = []
			    isatabfile.assays_factors[k]['accession'] =  []
			isatabfile.assays_factors[k]['value'].append(chemical)
			isatabfile.assays_factors[k]['ref'].append(term.split(':')[0])
			isatabfile.assays_factors[k]['concentration'].append(concentration)
			isatabfile.assays_factors[k]['unit'].append(units)
			isatabfile.assays_factors[k]['accession'].append(term)
			chemicals.add(k.replace('-',' '))
			chemicals.add(k.replace('-',' ') + ' concentration')
		    else:
			if k in chemicals:
			    isatabfile.assays_factors[k]['value'].append('')
			    isatabfile.assays_factors[k]['ref'].append('')
			    isatabfile.assays_factors[k]['concentration'].append('')
			    isatabfile.assays_factors[k]['unit'].append('')
			    isatabfile.assays_factors[k]['accession'].append('')
	        for p in general_factors:
		    if p[0] in isatabfile.assays_dict:
			isatabfile.assays_dict[p[0]].append(p[1])
		    else:
			isatabfile.assays_dict[p[0]] = [p[1]]
	    for p in protocols:
		if p in ontology.protocols:
		    term = ontology.protocols[p]
		else:
		    term = p
		isatabfile.investigation_dict['Study Protocol Name'].append(p.replace('-',' '))
		isatabfile.investigation_dict['Study Protocol Type Term Source REF'].append('OBI')
		isatabfile.investigation_dict['Study Assay Measurement Type'].append(p.replace('-',' '))
		isatabfile.investigation_dict['Study Assay Measurement Type Term Accession Number'].append(term)
		isatabfile.investigation_dict['Study Assay Measurement Type Term Source REF'].append(term.split(':')[0])
		isatabfile.investigation_dict['Study Assay Technology Type'].append(tech)
		isatabfile.investigation_dict['Study Assay File Name'].append('study-assay.txt')
	    for ch in chemicals:
		isatabfile.investigation_dict['Study Factor Name'].append(ch)
		isatabfile.investigation_dict['Study Factor Type'].append('Compound')
		isatabfile.investigation_dict['Study Factor Type Term Source REF'].append('CHEBI')
	    return isatabfile

		    
"""
For the ISATAB format for chemical footprinting experiments
"""


investigation_keys = ['ONTOLOGY SOURCE REFERENCE','Term Source Name','Term Source File','Term Source Version','Term Source Description','INVESTIGATION','Investigation Identifier','Investigation Title','Investigation Description','Investigation Submission Date','Investigation Public Release Date','INVESTIGATION PUBLICATIONS','Investigation PubMed ID','Investigation Publication DOI','Investigation Publication Author list','Investigation Publication Title','Investigation Publication Status','Investigation Publication Status Term Accession Number','Investigation Publication Status Term Source REF','INVESTIGATION CONTACTS','Investigation Person Last Name','Investigation Person First Name','Investigation Person Mid Initials','Investigation Person Email','Investigation Person Phone','Investigation Person Fax','Investigation Person Address','Investigation Person Affiliation','Investigation Person Roles','Investigation Person Roles Term Accession Number','Investigation Person Roles Term Source REF','STUDY','Study Identifier','Study Title','Study Submission Date','Study Public Release Date','Study Description','Study File Name','STUDY DESIGN DESCRIPTORS','Study Design Type','Study Design Type Term Accession Number','Study Design Type Term Source REF','STUDY PUBLICATIONS','Study PubMed ID','Study Publication DOI','Study Publication Author list','Study Publication Title','Study Publication Status','Study Publication Status Term Accession Number','Study Publication Status Term Source REF','STUDY FACTORS','Study Factor Name','Study Factor Type','Study Factor Type Term Accession Number','Study Factor Type Term Source REF','Study Assay Measurement Type','Study Assay Measurement Type Term Accession Number','Study Assay Measurement Type Term Source REF','Study Assay Technology Type','Study Assay Technology Type Term Accession Number','Study Assay Technology Type Term Source REF','Study Assay Technology Platform','Study Assay File Name','STUDY PROTOCOLS','Study Protocol Name','Study Protocol Type','Study Protocol Type Term Accession Number','Study Protocol Type Term Source REF','Study Protocol Description','Study Protocol URI','Study Protocol Version','Study Protocol Parameters Name','Study Protocol Parameters Name Term Accession Number','Study Protocol Parameters Name Term Source REF','Study Protocol Components Name','Study Protocol Components Type','Study Protocol Components Type Term Accession Number','Study Protocol Components Type Term Source REF','STUDY CONTACTS','Study Person Last Name','Study Person First Name','Study Person Mid Initials','Study Person Email','Study Person Phone','Study Person Fax','Study Person Address','Study Person Affiliation','Study Person Roles', 'Study Person Roles Term Accession Number', 'Study Person Roles Term Source REF']
assays_optional_columns = ['Factor Value[chemical]', 'Term Source REF[chemical]', 'Term Accession Number[chemical]', 'Factor Value[chemical concentration]', 'Unit[chemical]']
class ISATABFile:
    def __init__(self):
	self.assays_keys = ['Source Name','Characteristics[Nucleotide Sequence]','Characteristics[Nucleotide Type]','Characteristics[RNA Production]','Protocol REF','Parameter Value[Data Starts at Sequence Position]', 'Term Source REF','Term Accession Number', 'Assay Name', 'Raw Data File','Performer','Date']
	self.investigation_dict = {}
	self.assays_dict = {}
	self.assays_dict['Source Name'] = []
	self.assays_dict['Characteristics[Nucleotide Sequence]'] = []
	self.assays_dict['Characteristics[Nucleotide Type]'] = []
	self.assays_dict['Characteristics[RNA Production]'] = []
	self.assays_dict['Protocol REF'] = []
	self.assays_dict['Parameter Value[Data Starts at Sequence Position]'] = []
	self.assays_dict['Assay Name'] = []
	self.assays_dict['Raw Data File'] = []
	self.assays_dict['Performer'] = []
	self.assays_dict['Date'] = []
	self.assays_dict['Term Source REF'] = []
	self.assays_dict['Term Accession Number'] = []
	self.investigation_dict['ONTOLOGY SOURCE REFERENCE'] = []
	self.investigation_dict['Term Source Name'] = ['OBI','CHEBI', 'UO']
	self.investigation_dict['Term Source File'] = [] 
	self.investigation_dict['Term Source Version'] = ['v 1.26', 'v 1.26', 'v 1.26']
	self.investigation_dict['Term Source Description'] = ['Ontology for Biomedical Investigations', 'Chemical Entity of Biological Interest', 'Unit Ontology']
	self.investigation_dict['INVESTIGATION'] = []
	self.investigation_dict['Investigation Identifier'] = []
	self.investigation_dict['Investigation Title'] = []
	self.investigation_dict['Investigation Description'] = []
	self.investigation_dict['Investigation Submission Date'] = []
	self.investigation_dict['Investigation Public Release Date'] = []
	self.investigation_dict['INVESTIGATION PUBLICATIONS'] = []
	self.investigation_dict['Investigation PubMed ID'] = []
	self.investigation_dict['Investigation Publication DOI'] = []
	self.investigation_dict['Investigation Publication Author list'] = []
	self.investigation_dict['Investigation Publication Title'] = []
	self.investigation_dict['Investigation Publication Status'] = []
	self.investigation_dict['Investigation Publication Status Term Accession Number'] = []
	self.investigation_dict['Investigation Publication Status Term Source REF'] = []
	self.investigation_dict['INVESTIGATION CONTACTS'] = []
	self.investigation_dict['Investigation Person Last Name'] = []
	self.investigation_dict['Investigation Person First Name'] = []
	self.investigation_dict['Investigation Person Mid Initials'] = []
	self.investigation_dict['Investigation Person Email'] = []
	self.investigation_dict['Investigation Person Phone'] = []
	self.investigation_dict['Investigation Person Fax'] = []
	self.investigation_dict['Investigation Person Address'] = []
	self.investigation_dict['Investigation Person Affiliation'] = []
	self.investigation_dict['Investigation Person Roles'] = []
	self.investigation_dict['Investigation Person Roles Term Accession Number'] = []
	self.investigation_dict['Investigation Person Roles Term Source REF'] = []
	self.investigation_dict['STUDY'] = []
	self.investigation_dict['Study Identifier'] = []
	self.investigation_dict['Study Title'] = []
	self.investigation_dict['Study Submission Date'] = []
	self.investigation_dict['Study Public Release Date'] = []
	self.investigation_dict['Study Description'] = []
	self.investigation_dict['Study File Name'] = []
	self.investigation_dict['STUDY DESIGN DESCRIPTORS'] = []
	self.investigation_dict['Study Design Type'] = ['SNRNASM']
	self.investigation_dict['Study Design Type Term Accession Number'] = ['OBI']
	self.investigation_dict['Study Design Type Term Source REF'] = []
	self.investigation_dict['STUDY PUBLICATIONS'] = []
	self.investigation_dict['Study PubMed ID'] = []
	self.investigation_dict['Study Publication DOI'] = []
	self.investigation_dict['Study Publication Author list'] = []
	self.investigation_dict['Study Publication Title'] = []
	self.investigation_dict['Study Publication Status'] = []
	self.investigation_dict['Study Publication Status Term Accession Number'] = []
	self.investigation_dict['Study Publication Status Term Source REF'] = []
	self.investigation_dict['STUDY FACTORS'] = []
	self.investigation_dict['Study Factor Name'] = []
	self.investigation_dict['Study Factor Type'] = []
	self.investigation_dict['Study Factor Type Term Accession Number'] = []
	self.investigation_dict['Study Factor Type Term Source REF'] = []
	self.investigation_dict['Study Assay Measurement Type'] = []
	self.investigation_dict['Study Assay Measurement Type Term Accession Number'] = []
	self.investigation_dict['Study Assay Measurement Type Term Source REF'] = []
	self.investigation_dict['Study Assay Technology Type'] = []
	self.investigation_dict['Study Assay Technology Type Term Accession Number'] = []
	self.investigation_dict['Study Assay Technology Type Term Source REF'] = []
	self.investigation_dict['Study Assay Technology Platform'] = []
	self.investigation_dict['Study Assay File Name'] = []
	self.investigation_dict['STUDY PROTOCOLS'] = []
	self.investigation_dict['Study Protocol Name'] = []
	self.investigation_dict['Study Protocol Type'] = []
	self.investigation_dict['Study Protocol Type Term Accession Number'] = []
	self.investigation_dict['Study Protocol Type Term Source REF'] = []
	self.investigation_dict['Study Protocol Description'] = []
	self.investigation_dict['Study Protocol URI'] = []
	self.investigation_dict['Study Protocol Version'] = []
	self.investigation_dict['Study Protocol Parameters Name'] = ['data starts at sequence position']
	self.investigation_dict['Study Protocol Parameters Name Term Accession Number'] = []
	self.investigation_dict['Study Protocol Parameters Name Term Source REF'] = []
	self.investigation_dict['Study Protocol Components Name'] = []
	self.investigation_dict['Study Protocol Components Type'] = []
	self.investigation_dict['Study Protocol Components Type Term Accession Number'] = []
	self.investigation_dict['Study Protocol Components Type Term Source REF'] = []
	self.investigation_dict['STUDY CONTACTS'] = []
	self.investigation_dict['Study Person Last Name'] = []
	self.investigation_dict['Study Person First Name'] = []
	self.investigation_dict['Study Person Mid Initials'] = []
	self.investigation_dict['Study Person Email'] = []
	self.investigation_dict['Study Person Phone'] = []
	self.investigation_dict['Study Person Fax'] = []
	self.investigation_dict['Study Person Address'] = []
	self.investigation_dict['Study Person Affiliation'] = []
	self.investigation_dict['Study Person Roles'] = []
	self.investigation_dict['Study Person Roles Term Accession Number'] = []
	self.investigation_dict['Study Person Roles Term Source REF'] = []
	self.sample_id_name_map = {}
	self.data_id_order = []
	self.assays_factors = {}
	self.data = {}

    def save(self, name, type='xls'):
        self.name = name
        global investigation_keys
	if type == 'dir':
	    if not os.path.exists(name):
		os.mkdir(name)
	    investigationfile = open(name+'/investigation.txt', 'w')
	    assaysfile = open(name+'/study-assay.txt', 'w')
	    datamatrixfile = open(name+'/datamatrix.txt', 'w')
	    for k in self.assays_keys:
		assaysfile.write(k + '\t')
	    for k in self.assays_factors:
		assaysfile.write('Factor Value[%s]\t' % k.replace('-',' '))
		assaysfile.write('Term Source REF\t')
		assaysfile.write('Term Accession Number\t')
		assaysfile.write('Factor Value[%s concentration]\t' % k.replace('-',' '))
		assaysfile.write('Unit\t')
	    assaysfile.write('\n')
	    for i in range(len(self.assays_dict.values()[0])):
		line = ''
		for k in self.assays_keys:
		    if k in self.assays_dict:
			if len(self.assays_dict[k]) <= i:
			    line += '\t'
			else:
			    line += self.assays_dict[k][i] + '\t' 
		for k in self.assays_factors:
		    line += self.assays_factors[k]['value'][i] + '\t'	
		    line += self.assays_factors[k]['ref'][i] + '\t'	
		    line += self.assays_factors[k]['accession'][i] + '\t'	
		    line += self.assays_factors[k]['concentration'][i] + '\t'	
		    line += self.assays_factors[k]['unit'][i] + '\t'	
		assaysfile.write(line.strip('\t') + '\n')
	    for k in investigation_keys:
		line = k + '\t'
		for i in range(len(self.investigation_dict[k])):
		    line += self.investigation_dict[k][i] + '\t' 
		investigationfile.write(line.strip('\t') + '\n')
	    maxlen = max([len(x) for x in self.data.values()])
	    datamatrixfile.write('\t'.join(self.assays_dict['Source Name']) + '\n')
	    for i in range(maxlen):
		datamatrixfile.write('\t'.join(['' if i >= len(self.data[k]) else str(self.data[k][i]) \
						for k in self.data_id_order]) + '\n')
	    assaysfile.close()
	    datamatrixfile.close()
	    investigationfile.close()
	elif type == 'xls':
	    assayrow = 0
	    invrow = 0
	    datarow = 0
	    wb =xlwt.Workbook()
	    investigationsh = wb.add_sheet('investigation')
	    assayssh = wb.add_sheet('study-assay')
	    datamatrixsh = wb.add_sheet('datamatrix')
	    for i, k in enumerate(self.assays_keys):
		assayssh.write(assayrow, i, k)
	    for k in self.assays_factors:
		i += 1
		assayssh.write(assayrow, i , 'Factor Value[%s]\t' % k.replace('-',' '))
		i += 1
		assayssh.write(assayrow, i, 'Term Source REF\t')
		i += 1
		assayssh.write(assayrow, i, 'Term Accession Number\t')
		i += 1
		assayssh.write(assayrow, i, 'Factor Value[%s concentration]\t' % k.replace('-',' '))
		i += 1
		assayssh.write(assayrow, i, 'Unit\t')
	    
	    assayrow += 1
	    for i in range(len(self.assays_dict.values()[0])):
		line = []
		for k in self.assays_keys:
		    if k in self.assays_dict:
			if len(self.assays_dict[k]) <= i:
			    line.append('')
			else:
			    line.append(self.assays_dict[k][i]) 
	        if i < len(self.assays_factors):
		    for k in self.assays_factors:
			line .append( self.assays_factors[k]['value'][i])	
			line .append( self.assays_factors[k]['ref'][i])	
			line .append( self.assays_factors[k]['accession'][i])	
			line .append( self.assays_factors[k]['concentration'][i])	
			line .append( self.assays_factors[k]['unit'][i])	
		for j in range(len(line)):
		    assayssh.write(assayrow, j, line[j])
		assayrow += 1
	    for i, k in enumerate(investigation_keys):
		line = [k]
		for j in range(len(self.investigation_dict[k])):
		    line.append(self.investigation_dict[k][j]) 
		for j in range(len(line)):
		    investigationsh.write(invrow, j, line[j])
		invrow += 1
	    maxlen = max([len(x) for x in self.data.values()])
	    for i, k in enumerate(self.assays_dict['Source Name']):
		datamatrixsh.write(datarow, i, k)
	    datarow += 1
	    for i in range(maxlen):

		for j, k in enumerate(['' if i >= len(self.data[k]) else str(self.data[k][i]) \
						for k in self.data_id_order]):
		    datamatrixsh.write(datarow, j, k)
		datarow += 1
            wb.save(name)

	else:
	    print 'Unrecognized type %s to save isatab file' % type

    def load(self, name, type='xls'):
	self.name = name
        if type == 'dir':
	    investigationfile = open(name+'/investigation.txt')
	    for l in investigationfile.readlines():
		fields = l.strip().split('\t')
		self.investigation_dict[fields[0]] = fields[1:]
	    assaysfile = open(name+'/'+self.investigation_dict['Study Assay File Name'][0])
	    assays_keys = assaysfile.readline().strip().split('\t')
	    l = assaysfile.readline()
	    while l:
		for i, f in enumerate(l.strip().split('\t')):
		    if assays_keys[i] in self.assays_dict:
			self.assays_dict[assays_keys[i]].append(f)
		    else:
			self.assays_dict[assays_keys[i]] = [f]
		l = assaysfile.readline()
	    datamatrixfile = open(name+'/'+self.assays_dict['Raw Data File'][0])
	    data_keys = datamatrixfile.readline().strip().split('\t')
	    for i, k in enumerate(data_keys):
		self.data[k + '_' + str(i+1)] = []
		self.sample_id_name_map[k + '_' + str(i+1)] = k
		self.data_id_order.append(k+ '_' + str(i+1))
	    l = datamatrixfile.readline()
	    while l:
		for i, f in enumerate(l.strip().split('\t')):
		    if f != '':
			self.data[data_keys[i] + '_' + str(i+1)].append(float(f))
		l = datamatrixfile.readline()
	    assaysfile.close()
	    datamatrixfile.close()
	    investigationfile.close()
	elif type == 'xls':
	    wb = xlrd.open_workbook(name)
	    investigationsh = wb.sheet_by_name('investigation')
	    for j in range(investigationsh.nrows):
		fields = investigationsh.row_values(j)
		self.investigation_dict[fields[0]] = fields[1:]
	    try:
		assayssh = wb.sheet_by_name(self.investigation_dict['Study Assay File Name'][0].replace('.txt',''))
	    except xlrd.biffh.XLRDError:
		assayssh = wb.sheet_by_name(self.investigation_dict['Study File Name'][0].replace('.txt',''))
	    assays_keys = assayssh.row_values(0)
	    for j in range(1, assayssh.nrows):
		l = assayssh.row_values(j)
		for i, f in enumerate(l):
		    if assays_keys[i] in self.assays_dict:
			self.assays_dict[assays_keys[i]].append(f)
		    else:
			self.assays_dict[assays_keys[i]] = [f]
	    datamatrixsh = wb.sheet_by_name(self.assays_dict['Raw Data File'][0].replace('.txt',''))
	    data_keys = datamatrixsh.row_values(0)
	    for i, k in enumerate(data_keys):
		self.data[k + '_' + str(i+1)] = []
		self.sample_id_name_map[k + '_' + str(i+1)] = k
		self.data_id_order.append(k+ '_' + str(i+1))
	    for j in range(1, datamatrixsh.nrows):
		l = datamatrixsh.row_values(j)
		for i, f in enumerate(l):
		    if f != '':
			self.data[data_keys[i] + '_' + str(i+1)].append(float(f))
	else:
	    print 'Unrecognized type %s for loading isatab file' % type

    def validate(self):
        messages = []
        def check_terms(d, prefix, ontdict):
	    m = []
	    for i, t in enumerate(d[prefix]):
	        if not t.replace(' ','-') in ontdict:
		    if t.strip() == '':
			pass
		    else:
			m.append( 'WARNING! For %s, term %s is unknown for its respective ontology' % (prefix, t))
		else:
		    r = d[prefix + ' Term Accession Number'][i].strip().replace('_',':')
		    if ontdict[t.replace(' ','-')] != r:
			m.append( 'WARNING! For %s, term %s and accession number %s do not match' % (prefix, t, r))
	    for i, t in enumerate(d[prefix + ' Term Accession Number']):
		if not d[prefix + ' Term Source REF'][i] in t:
		    m.append( 'WARNING! For %s, accession number and REF ontology do not match.' % prefix)
	    return m

        if len(self.data) != len(self.assays_dict['Source Name']):
	    messages.append( 'WARNING! Number of samples in assays and data file do not match')
        for k in self.data:
	    if not self.sample_id_name_map[k] in self.assays_dict['Source Name']:
		messages.append( 'WARNING! Sample %s in data file not referenced from assays file' % k)
	for i, k in enumerate(self.assays_dict['Source Name']):
	    if not k + '_' + str(i) in self.data:
		messages.append( 'WARNING! Sample %s in assays file is missing from data file' % k)
        for p in self.investigation_dict['Study Protocol Name']:
	    if p.strip() != '' and not p.strip().lower() in [x.strip().lower() for x in self.assays_dict['Protocol REF']]:
	        messages.append( 'WARNING! Protocol %s in investigation file is not in assays file' % p)
        # Checks for ontology term consistency
        #messages = messages + check_terms(self.investigation_dict, 'Study Assay Measurement Type', ontology.protocols) 
        #messages = messages + check_terms(self.investigation_dict, 'Study Factor Type', ontology.chemicals) 
	for i, seq in enumerate(self.assays_dict['Characteristics[Nucleotide Sequence]']):
	    if i >= len(self.assays_dict['Source Name']):
		messages.append( 'ERROR! Cannot continue validation as list of source names and sequences do not match in number!')
		return messages
	    for k in ['Source Name', 'Parameter Value[Data Starts at Sequence Position]', 'Characteristics[Nucleotide Type]']:
		stop = False
		if len(self.assays_dict[k]) == 0:
		    messages.append('ERROR! Cannot continue validation as "%s" column is not present in assays tab! (check any spelling inconsistencies in column names)' % k)
		    stop = True
		if stop:
		    return messages
 
            sn = self.assays_dict['Source Name'][i]
            idn = sn + '_' + str(i)
	    m = int(self.assays_dict['Parameter Value[Data Starts at Sequence Position]'][i])
	    if idn in self.data and len(self.data[idn]) != len(seq) - m + 1:
		    messages.append('WARNING! Length of data lane %s [%s] and length of respective sequence [%s] do not match' % (sn, len(self.data[idn]), len(seq) - m + 1))
	    chartype = self.assays_dict['Characteristics[Nucleotide Type]'][i]
	    if chartype == 'RNA' and 'T' in seq:
		messages.append( 'WARNING! Sequence for %s specified as RNA but looks like DNA.' % n)
	    if chartype == 'DNA' and 'U' in seq:
		messages.append( 'WARNING! Sequence for %s specified as DNA but looks like RNA.' % n)
        return messages

    def toRDAT(self):
	rdatfile = RDATFile()
	rdatfile.comments = ''
	general_annotations = defaultdict(list)
	for i, name in enumerate(self.assays_dict['Assay Name']):
	    rdatfile.constructs[name] = RDATSection()
	    d = RDATSection()
	    c = rdatfile.constructs[name]
	    c.name = name
	    c.annotations = {}
	    for k, v in general_annotations.iteritems():
		c.annotations[k] = v
	    d.values = self.data[name]
	    rdatfile.values[name] = [d.values]
	    d.annotations = {}
	    d.xsel = []
	    d.errors = []
	    d.trace = []
	    d.seqpos = []
	    c.sequence = self.assays_dict['Characteristics[Nucleotide Sequence]'][i]
	    c.seqpos = range(len(c.sequence))
	    c.xsel = []
	    c.values = []
	    c.traces = []
	    c.mutpos = []
	    c.data_types = []
	    c.data = [d]
	    c.offset = 0
	    c.structure = '.'*len(c.sequence)
	rdatfile.loaded = True
	return rdatfile

