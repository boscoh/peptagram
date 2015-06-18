DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
cd $DIR

python do_xtandem_peptagram.py $*

if [ "$(uname)" == "Darwin" ]; then
  echo -n -e "]0;mac_xtandem_peptagram.command.shell"
  osascript -e 'tell application "Terminal" to close (every window whose name contains "mac_xtandem_peptagram.command.shell")' &
fi
