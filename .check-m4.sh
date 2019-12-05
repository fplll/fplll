#!/bin/bash

for f in m4/*.m4 ; do
	oldserial=$(grep "^#serial" $f | grep -o "[0-9]*" | head -n1)
	if [ -z $oldserial ]; then
		continue
	fi

	echo "Checking for newer version of: $f"
	if wget -t 1 -T 0.5 -O $f.tmp "https://raw.githubusercontent.com/autoconf-archive/autoconf-archive/master/$f" >/dev/null 2>/dev/null ; then
		newserial=$(grep "^#serial" $f.tmp | grep -o "[0-9]*" | head -n1)
		rm $f.tmp
	else
		echo "WGET error"
		rm $f.tmp
		continue
	fi
	if [ -z $newserial ]; then
		echo "Error checking newer version"
		continue
	fi

	if [ $oldserial -eq $newserial ]; then
		echo "Current version is latest version"
	else
		echo "*** UPDATE AVAILABLE for $f: version $oldserial (current) => $newserial (latest) ***"
	fi
done
