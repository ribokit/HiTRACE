from setuptools import setup

setup(name='RDATkit',
        version='0.5',
        description='RNA dataset toolkit',
        author='Pablo Cordero',
        author_email='tsuname@stanford.edu',
        url='http://rmdb.stanford.edu/rdatkit',
        packages=['rdatkit', 'rdatkit.likelihood', 'rdatkit.mutate_and_map'],
        install_requires=['xlrd', 'xlwt']
        )
