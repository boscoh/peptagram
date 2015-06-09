DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
cd $DIR
python pilot_peptagram.py $*
if [ "$(uname)" == "Darwin" ]; then
  echo -n -e "]0;pilot_peptagram.command"
  osascript -e 'tell application "Terminal" to close (every window whose name contains "pilot_peptagram.command")' &
fi