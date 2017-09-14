#!/bin/bash

function is_file() {
    if [ ! -f $1 ]; then
	echo "File not found: $1"
	return 1
    fi
    return 0
}

function check_scan() {
    BOUT=$1
    is_file $BOUT || return
    if grep -q -E "[0-9]+ (bug|bugs) found" $BOUT; then
       DIR=`tail -n 1 $BOUT | egrep -o "[0-9]+-[0-9]+-[0-9]+-[0-9]+-[0-9]+-[0-9]"`;
       echo; echo "View report at http://`hostname -f`/scan-builds/$DIR";
       return 2;
    fi
}

function grep_log() {
    while read LOG; do

    	is_file $LOG || continue

        if ! grep '^# FAIL: *0$' $LOG &> /dev/null
        then
	  echo "===== $LOG ====="
	  tail -n 200 $LOG
	  exit 2
        fi
     done
}

function check_check() {
    export -f grep_log
    export -f is_file
    find . -name "test-suite.log" | grep_log
}

