DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
cd $DIR
python mascot_peptagram.py $*
if [ "$(uname)" == "Darwin" ]; then
  echo -n -e "]0;mascot_peptagram"
  osascript -e 'tell application "Terminal" to close (every window whose name contains "mascot_peptagram")' &
fi
