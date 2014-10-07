#!/usr/bin/python
# -*- coding: utf-8 -*-
from optparse import OptionParser

def do_plot(filename):
    import matplotlib.pyplot as plt
    x1,x2,y1,y2 = None,None,None,None
    for i,line in enumerate(open(filename).readlines()):

        plt.clf()
        line = line.strip()
        name, data = line.split(":")
        data = [d for d in data.split(" ") if d]
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