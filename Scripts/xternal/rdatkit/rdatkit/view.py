import settings
import pdb
import os
from secondary_structure import SecondaryStructure

class VARNA:
    
    @classmethod
    def get_colorMapStyle(self, values, sequence=''):
        if len(sequence) > len(values):
	    return '-0.001:#C0C0C0,0:#FFFFFF;0.1:#FFFFFF,0.8:#FF8800;1:#FF0000'
        if reduce(lambda x,y: x or y, [x < 0 for x in values]):
	    return '-0.001:#C0C0C0,0:#FFFFFF;0.1:#FFFFFF,0.8:#FF8800;1:#FF0000'
	else:
	    return '0:#FFFFFF;0.1:#FFFFFF,0.8:#FF8800;1:#FF0000'
    
    @classmethod
    def cmd(self, sequence, structure, outfile, options={}):
        option_str = ''
	for key in options:
	    val = options[key]
	    if type(val) == list:
		argval = str(val).strip('[]').replace('L', '').replace('u','')
	    else:
		argval = str(val)
	    option_str += '-%s "%s" ' % (key, argval)
	os.popen('java -cp %s  fr.orsay.lri.varna.applications.VARNAcmd -sequenceDBN %s -structureDBN "%s" %s -o %s' %\
		  (settings.VARNA, sequence, structure, option_str, outfile))

    def __init__(self, sequences=[], structures=[], mapping_data=[]):
        self.sequences = sequences
	self.structures = structures
        self.mapping_data = mapping_data
	self.rows = 1
	self.columns = 1
	self.width = 522
	self.height = 200
	self.annotation_font_size = 9
    
    def get_values(self, att):
        if att == 'structures':
	    return [struct.dbn for struct in self.structures]
        return self.__dict__[att]
    
    def get_frames(self, overlap_structures):
        res = -1
        for key, val in self.__dict__.items():
	    if overlap_structures and key == 'structures':
		continue
	    if type(val) == list:
		res = max(res, len(val))
	return res

    def render(self, base_annotations=[], annotation_by_helix=False, helix_function=(lambda x,y:x), map_data_function=(lambda x:x), overlap_structures=False, reference_structure=SecondaryStructure(), annotation_def_val=''):
        struct_string = ''
	applet_string = '<applet  code="VARNA.class" codebase="http://rmdb.stanford.edu/site_media/bin"\n'
	applet_string += 'archive="VARNA.jar" width="%d" height="%d">\n' % (self.width, self.height)
	frames = self.get_frames(overlap_structures=overlap_structures)
	if overlap_structures:
	    struct_string += '<param name="structureDBN" value="%s"/>\n' % self.structures[0].dbn
	    bps = self.structures[0].base_pairs()
	    struct_string += '<param name="auxBPs" value="'
	    for i in range(1, len(self.structures)):
		for bp in self.structures[i].base_pairs():
		    if bp not in bps:
			bps.append(bp)
			struct_string += '(%s,%s):edge5=s,edge3=h,stericity=cis;' % (bp[0]+1, bp[1]+1)
	    struct_string += '"/>\n'
	base_annotation_string = ''
	for i, ba in enumerate(base_annotations):
	    if len(base_annotations) > 1:
		idx = str(i+1)
	    else:
		idx = ''
	    base_annotation_string += '<param name="annotations%s" value=\n\t"' % idx
	    if annotation_by_helix:
		for helix in self.structures[i].helices():
		    anchor = helix[0][0] + 1 + (helix[-1][0] - helix[0][0])/2 
		    annotation_value = ba.values()[0]
		    for bp in helix:
			if bp in ba:
			    nextval = ba[bp]
			    annotation_value = helix_function(annotation_value, nextval)
		    base_annotation_string += '%s:type=L,anchor=%d,size=%d;\n' % (annotation_value, anchor, self.annotation_font_size)
	    else:
		for bp in base_annotations:
		    base_annotation_string += '%s:type=B,anchor=%d,size=%d;\n' % (base_annotations[bp], bp[0], self.annotation_font_size)	
            base_annotation_string = base_annotation_string.strip() + '"/>\n'
	    applet_string += base_annotation_string
	if self.rows * self.columns != frames:
	    rows = frames
	    columns = 1
	else:
	    rows = self.rows
	    columns = self.columns
	param_string = ''
	param_string += '<param name="rows" value="%d"/>\n' % rows
	param_string += '<param name="columns" value="%d"/>\n' % columns
	do_default_colormapstyle = 'colorMapStyle' not in self.__dict__
	for att in self.__dict__:
	    if overlap_structures and att == 'structures':
		param_string += struct_string
	    if type(self.get_values(att)) == list:
		if att == 'structures':
		    if overlap_structures:
			continue
		    name = 'structureDBN'
		    if len(reference_structure) > 0:
			for i in range(len(self.structures)):
			    if len(self.structures) > 1:
				frame_idx = str(i+1)
			    else:
				frame_idx = ''
			    param_string += '<param name="%s%s" value="%s" />\n' % (name, frame_idx, self.structures[i].dbn)
			    auxbps_string = '<param name="auxBPs%s" value="' % (frame_idx)
			    refbps = reference_structure.base_pairs()
			    bps = self.structures[i].base_pairs()
			    for bp in bps:
				if bp not in refbps:
				    auxbps_string += '(%s,%s):color=#00CCFF;' % (bp[0]+1, bp[1]+1)
			    for bp in refbps:
				if bp not in bps:
				    auxbps_string += '(%s,%s):color=#FF6600;' % (bp[0]+1, bp[1]+1)
			    auxbps_string += '"/>\n'
			    param_string += auxbps_string
			continue
		elif att == 'sequences':
		    name = 'sequenceDBN'
		elif att == 'mapping_data':
		    name = 'colorMap'
		else:
		    name = att
		for i, val in enumerate(self.get_values(att)):
		    if frames > 1:
			frame_idx = str(i+1)
		    else:
			frame_idx = ''
                    if name == 'colorMap': 
			if len(self.sequences) == len(self.mapping_data):
			    idx = i
			else:
			    idx = 0
			if do_default_colormapstyle:
			    param_string += '<param name="colorMapStyle%s" value="%s" />\n' % (frame_idx, VARNA.get_colorMapStyle(val, sequence=self.sequences[idx]))
			val = [map_data_function(val[x]) if val[x] != None and val[x] > 0 else -0.001 for x in range(len(self.sequences[idx]))]
			val = str(val).strip('[]').replace(' ','')
		    param_string += '<param name="%s%s" value="%s" />\n' % (name, frame_idx, val)
	applet_string += param_string
	applet_string += '</applet>'
	return applet_string

	
	



