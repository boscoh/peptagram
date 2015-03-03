DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
cd $DIR
python prophet_peptagram.py $*
if [ "$(uname)" == "Darwin" ]; then
  echo -n -e "]0;prophet_peptagram.command"
  osascript -e 'tell application "Terminal" to close (every window whose name contains "prophet_peptagram.command")' &
fi
