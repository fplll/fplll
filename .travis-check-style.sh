#!/bin/bash

# Add LLVM toolchain
sudo add-apt-repository 'deb http://apt.llvm.org/trusty/ llvm-toolchain-trusty-3.9 main'
wget -O - http://llvm.org/apt/llvm-snapshot.gpg.key | sudo apt-key add -

# Install clang-format
sudo apt-get update -qq 
sudo apt-get install -qq -y clang-format-3.9

# If clang-format modifies files, the code does not adhere to the project standards
if [[ $(make codingstyle && git status -s) ]]; 
then 
	tput setaf 1;
	echo "Code does not adhere to the project standards. Run \"make codingstyle\".";
	exit 1;
else 
	tput setaf 2;
	echo "Code adheres to the project standards.";
	exit 0;
fi;
