DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
cd $DIR
python qc_peptagram.py $*
if [ "$(uname)" == "Darwin" ]; then
  echo -n -e "]0;qc_peptagram.command"
  osascript -e 'tell application "Terminal" to close (every window whose name contains "qc_peptagram.command")' &
fi
