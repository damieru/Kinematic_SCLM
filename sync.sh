#!/usr/bin/env bash

SCRIPT_DIR=$( cd "$( dirname "$0" )" && pwd )
readonly SCRIPT_DIR

function main()
{
  while true; do
    find "$SCRIPT_DIR" -type f |
        entr -d bash "$SCRIPT_DIR/entr_sync.sh" "$SCRIPT_DIR" 

    # This is needed to stop the script on Ctrl + C
    if [[ $? -eq 0 ]]; then
        exit
    fi
  done
}

main "$@"