#!/usr/bin/python
# -*- coding: utf-8 -*-
from optparse import OptionParser
import json

def do_plot(filename):
    import matplotlib.pyplot as plt
    x1,x2,y1,y2 = None,None,None,None
    i = 0
    filer = open(filename, "r")
    data_json = json.load(filer)
    filer.close()

    for d in data_json:

        plt.clf()
        name = d.get("step")
        if int(d.get("loop")) != -1:
            name += " {0}".format(d.get("loop"))
        if d.get("step") != "Input":
            name += " ({0} s)".format(d.get("time"))
        data = d.get("norms")
        data = map(float, data)
        
        plt.plot(range(len(data)), data)
        plt.xlabel(name)

        if i == 0:
            x1,x2,y1,y2 = plt.axis()
        else:
            axes = plt.gca()
            axes.set_xlim([x1,x2])
            axes.set_ylim([y1,y2])
        from os.path import basename
        name = basename(filename)
        plt.savefig("%s-%04d.png"%(name,i))
        i += 1

def main():
    parser = OptionParser(usage="""\
""")
    opts, args = parser.parse_args()

    try:
        filename = args[0]
    except IndexError:
        parser.print_help()
        return

    do_plot(filename)

if __name__ == '__main__':
    main()
