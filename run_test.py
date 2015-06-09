import os, glob
for py in glob.glob('*_peptagram.py'):
    if 'reorder' not in py:
        os.system('python {} test'.format(py))