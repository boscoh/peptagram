import glob
import os
import shutil

batch_template = """\
python %s %%1 %%2 %%3
"""

shell_template = """\
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
cd $DIR
python %s $*
"""

[os.remove(s) for s in glob.glob('*_peptagram')]
[os.remove(s) for s in glob.glob('*_peptagram.cmd')]

for py_script in glob.glob('*_peptagram.py'):
    batch = py_script.replace('.py', '.bat')
    open(batch, 'w').write(batch_template % py_script)
    os.system('chmod -x ' + batch)
    print batch

    shell = py_script.replace('.py', '.command')
    open(shell, 'w').write(shell_template % py_script)
    os.system('chmod +x ' + shell)
    print shell
