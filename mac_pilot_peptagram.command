DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
cd $DIR

python do_pilot_peptagram.py $*

if [ "$(uname)" == "Darwin" ]; then
  echo -n -e "]0;mac_pilot_peptagram.command.shell"
  osascript -e 'tell application "Terminal" to close (every window whose name contains "mac_pilot_peptagram.command.shell")' &
fi
