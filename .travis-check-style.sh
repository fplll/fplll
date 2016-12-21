#!/bin/bash

make check-style
if [[ $(git status -s) ]]; 
then
	git diff
	tput setaf 1;
	echo "Code does not adhere to the project standards. Run \"make check-style\".";
	exit 1;
else 
	tput setaf 2;
	echo "Code adheres to the project standards.";
	exit 0;
fi;
