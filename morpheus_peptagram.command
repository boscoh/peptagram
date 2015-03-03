DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
cd $DIR
python morpheus_peptagram.py $*
if [ "$(uname)" == "Darwin" ]; then
  echo -n -e "]0;morpheus_peptagram.command"
  osascript -e 'tell application "Terminal" to close (every window whose name contains "morpheus_peptagram.command")' &
fi
