DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
cd $DIR

python xtandem_peptagram.py $*

if [ "$(uname)" == "Darwin" ]; then
  echo -n -e "]0;xtandem_peptagram.command.shell"
  osascript -e 'tell application "Terminal" to close (every window whose name contains "xtandem_peptagram.command.shell")' &
fi
