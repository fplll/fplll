#!/usr/bin/perl
###########################################################################
# Copyright (C) 2013 Damien Stehle.
#
# This file is part of fplll. fplll is free software: you
# can redistribute it and/or modify it under the terms of the GNU Lesser
# General Public License as published by the Free Software Foundation,
# either version 2.1 of the License, or (at your option) any later version.
#
# fplll is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with fplll. If not, see <http://www.gnu.org/licenses/>. */
##########################################################################

$r=@ARGV[0];
$c=@ARGV[1];



$i=0;
while(<stdin>)
{
    while (s/^[^0123456789\-\+]*([\-]?[0-9]+)([^0-9].*)$/$2/)#$1=first number in chain, $2=queue of string -> $_. 
    {
	@v[$i++]=$1;
    }
}

die "bizarre : $i != $r * $c" if ($i!=$r*$c);
print "B:=RMatrixSpace(Integers(), $r, $c) ! [";
$k=0;
for ($i=0;$i<$r;$i++)
{
    for ($j=0;$j<$c;$j++)
    {
	print @v[$k++];
	print "," if ($i!=$r-1 || $j!=$c-1);
    }
}
print "];";
