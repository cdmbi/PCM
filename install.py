import subprocess
import os

try:
    import rdkit
except ImportError:
    print >> sys.stderr, 'rdkit not installed...'
    print >> sys.stderr, 'build rdkit...'
    # subprocess.call(['git', 'clone', 'https://github.com/rdkit/conda-rdkit.git'])

os.chdir('conda-rdkit')
subprocess.call(['conda', 'build', 'boost'])
subprocess.call(['conda', 'build', 'rdkit'])
