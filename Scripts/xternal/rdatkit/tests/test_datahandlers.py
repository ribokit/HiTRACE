import datahandlers
import sys

print 'Creating RDATFile object'
rdat = datahandlers.RDATFile()
print 'Loading RDAT file'
rdat.load(open(sys.argv[1]))
print 'Validating file'
rdat.validate()
print 'Converting to isatab'
isatab = rdat.toISATAB()
