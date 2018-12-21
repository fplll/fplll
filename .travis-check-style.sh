#!/bin/bash

if [ "$CI" = "true" ] && [ "$TRAVIS_OS_NAME" = "osx" ]; then
   echo "Skipping style check on OSX due to being unreliable.";
   exit 0;
fi

if [ "$CI" = "true" ] && [ "$COMPILER" = "g++-4.8" ]; then
   echo "Skipping style check when testing g++4.8 due to old clang-format-3.9.";
   exit 0;
fi


make check-style
if [[ $(git status -s) ]]; 
then
	git --no-pager diff
	tput setaf 1;
	echo "Code does not adhere to the project standards. Run \"make check-style\".";
	exit 1;
else 
	tput setaf 2;
	echo "Code adheres to the project standards.";
	exit 0;
fi;
