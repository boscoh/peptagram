DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
cd $DIR
python reorder_peptagram.py $*
if [ "$(uname)" == "Darwin" ]; then
  echo -n -e "]0;reorder_peptagram.command"
  osascript -e 'tell application "Terminal" to close (every window whose name contains "reorder_peptagram.command")' &
fi
