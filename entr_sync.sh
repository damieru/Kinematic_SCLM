#!/usr/bin/env bash

# Note: this file is intended to be executed by entr

readonly SYNC_FOLDER="$1"

function main()
{
    reset

    run_sync

    echo
    print_line
    date "+Date: %A, %x | Time: %T %Z"
}

function run_sync()
{
    rsync -avu --progress \
      --filter="- mfiles"     \
      --filter="- Plots"      \
      --filter="- Results"    \
      --filter="- .git"       \
      --filter="- .gitignore" \
      --filter="- .vscode"    \
      --filter="- *.mod"      \
      --filter="- *.o"        \
      --filter="- *.out"      \
      "$SYNC_FOLDER" funk:/STORAGE/DATA/dalbuquerque/

}


function print_line()
{
  echo "=================================================="
}

main "$@"