from datahandlers import *
import sys

rdat = RDATFile()
rdat.load(open(sys.argv[1]))
isatab = rdat.toISATAB()
isatab.save(sys.argv[1][:sys.argv[1].find('.')] + '.xls', type='xls')
