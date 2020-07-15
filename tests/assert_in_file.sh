#!/bin/bash

if [ $# -ne 3 ]; then
    echo $0: usage: assert_in_file.sh filepath query inverse
    exit 1
fi

if [[ "$3" == "1" ]]; then
	echo -e "\nTesting if '$2' is not in '$1'"
	if grep -q "$2" "$1"; then
		echo -e "\ERROR: '$2' is in '$1'\nContents:\n"
	    cat "$1"
	    exit 1
	else
	    echo -e "\nSUCCESS: '$2' not in '$1'\n"
	    exit 0
	fi
else
	echo -e "\nTesting if '$2' is in '$1'"
	if grep -q "$2" "$1"; then
	    echo -e "\nSUCCESS: '$2' is in '$1'\n"
	    exit 0
	else
	    echo -e "\nERROR: '$2' not in '$1'\nContents:\n"
	    cat "$1"
	    exit 1
	fi
fi
