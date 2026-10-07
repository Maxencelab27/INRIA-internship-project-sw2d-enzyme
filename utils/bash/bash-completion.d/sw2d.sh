#!/usr/bin/env bash

# BASH WARNING
# This autocompletion routine is made to run under bash version >=4.
# Please make sure your are using a comptaible version.
# See https://www.shell-tips.com/mac/upgrade-bash/ for MacOS users.
# This file should be copied in your autocompletion directory, e.g. /opt/homebrew/etc/bash_completion.d

_sw2d() {
  COMPREPLY=()

  # All possible first values in command line
  local SERVICES=("-b" "--build" "-d" "--doc" "-h" "--help" "-g" "--graphics" "-r" "--run" "-t" "--test")

  # declare an associative array for options (values should be read from outside)
  declare -A ACTIONS
  ACTIONS["-b"]="clean Coverage Debug Minisizerel Release RelWithDebInfo"
  ACTIONS["--build"]=${ACTIONS["-b"]}
  ACTIONS["-d"]="all clean html pdf"
  ACTIONS["--doc"]=${ACTIONS["-d"]}
  ACTIONS["-h"]=""
  ACTIONS["--help"]=${ACTIONS["-h"]}
   ACTIONS["-g"]=""
  ACTIONS["--graphics"]=${ACTIONS["-g"]}
  ACTIONS["-r"]="converter porositymapper postprocess solver zone"
  ACTIONS["--run"]=${ACTIONS["-r"]}
  ACTIONS["-t"]="clean converter ddp dip full liswes porositymapper postprocess solver sp swes unit zone"
  ACTIONS["--test"]=${ACTIONS["-r"]}


  # All possible SERVICES and ACTIONS
  KEYWORDS="sw2d -b --build clean Debug Minisizerel Release RelWithDebInfo -d --doc html -h --help -g --graphics -r --run -t --test"

  # current word being autocompleted
  local cur=${COMP_WORDS[COMP_CWORD]}

  if [ ${ACTIONS[$3]+1} ] ; then
    # If previous arg is a key of ACTIONS (so it is a service).
    # It means that we must display action choices
    compopt +o nospace
    COMPREPLY=( `compgen -W "${ACTIONS[$3]}" -- $cur` )
  elif [[ "${ACTIONS[-b]}" == *"$3"* ]] || [[ "${ACTIONS[-d]}" == *"$3"* ]] || [[ "${ACTIONS[-r]}" == *"$3"* ]] || [[ "${ACTIONS[-t]}" == *"$3"* ]] ; then
    # If previous arg is one of the actions of -b, -d, -r or -t or previous arg is a final option
    # We are at the end of the command and nothing should be available
    compopt +o nospace
    COMPREPLY=()
  elif [[ ${COMP_WORDS[COMP_CWORD]::1} == "-"  ]] ; then
    # We are currently writing a service
    compopt +o nospace
    COMPREPLY=(`compgen -W "${SERVICES[*]}" -- $cur`)
  elif [[ $3 == sw2d ]] ; then
    # Previous arg is the function name, sw2d
    compopt +o nospace
    COMPREPLY=(`compgen -W "${SERVICES[*]}" -d -- $cur`)
  fi
}
complete -F _sw2d -o nospace sw2d