import glob
import os
import shutil

batch_template = """\
start python %s -i
"""

shell_template = """\
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
cd $DIR

python %(py_script)s $*

if [ "$(uname)" == "Darwin" ]; then
  echo -n -e "\033]0;%(window_name)s\007"
  osascript -e 'tell application "Terminal" to close (every window whose name contains "%(window_name)s")' &
fi
"""

[os.remove(s) for s in glob.glob('*_peptagram')]
[os.remove(s) for s in glob.glob('*.command')]
[os.remove(s) for s in glob.glob('*.bat')]

for py_script in glob.glob('*_peptagram.py'):
    name = py_script.replace('.py', '').replace('do_', '')

    batch = 'win_' + name + '.bat'
    open(batch, 'w').write(batch_template % py_script)
    print batch

    shell = 'mac_' + name + '.command'
    sub = { 'py_script': py_script, 'window_name': shell+'.shell' }
    open(shell, 'w').write(shell_template % sub)
    os.system('chmod +x ' + shell)
    print shell
