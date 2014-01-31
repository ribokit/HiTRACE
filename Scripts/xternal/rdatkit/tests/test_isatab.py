import datahandlers
import sys
print 'Creating ISATABFile class'
file = datahandlers.ISATABFile()
print 'Loading files from directory'
if '.xls' in sys.argv[1]:
    type = 'xls'
else:
    type = 'dir'
file.load(sys.argv[1], type='xls')
print 'Validating'
file.validate()
print 'Writing'
file.save('test.xls', type='xls')
