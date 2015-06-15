# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to plot CG tail order parameters from CG and AA simulations
"""

import argparse
import sys

import matplotlib.pylab as plt
import numpy as np

import colors

def _plot_order(order_list,labels,figure):
    """
    Plot the distribution of a series in subplots
    """
    nchains = len(order_list[0])
    ncols = 2.0
    nrows = np.ceil(nchains / ncols)

    for ichain in range(nchains):
        a = figure.add_subplot(nrows, ncols, ichain + 1)
        for i,(simorder,label) in enumerate(zip(order_list,labels)):
            order = simorder[ichain].mean(axis=0)
            a.plot(np.arange(1,len(order)+1,1),order,'-*',label=label,color=colors.color(i))
        a.set_xticks(np.arange(1,len(order)+1,1))
        # Add legend and ticks
        if ichain == 0:
            a.legend(loc=1, fontsize=8)

def _read_series(filename, chains_count):
    """
    Read series of data output by md_order.py
    """
    data = []
    with open(filename, "r") as f:
        line = f.readline()
        while line:
            data.append(line.strip().split())
            line = f.readline()
    data = np.array(data, dtype=float)
    # Split the data for the different chains
    chainorder = []
    for i,cnt in enumerate(chains_count) :
        first = 0 if i == 0 else chains_count[-1]-1
        last = first + cnt 
        chainorder.append(data[:,first+1:last+1])
    return chainorder

if __name__ == '__main__':

    print " ".join(sys.argv)

    # Command-line input
    parser = argparse.ArgumentParser(description="Analyse CG order parameters")
    parser.add_argument('-f', '--file', nargs="+", help="the trajectory of the order parameters")
    parser.add_argument('-l', '--labels', nargs="+", help="the labels of the different trajectories")
    parser.add_argument('-c', '--chains',nargs="+", help="the carbon chains")
    parser.add_argument('-o', '--out', help="the output filename", default="order.png")
    args = parser.parse_args()

    # Parse the chain definition
    chains_beads = []
    chains_count = []
    for chain in args.chains :
        chains_beads.append(chain.split("-"))
        chains_count.append(len(chains_beads[-1])-1)

    # Read data, assume it was created by md_order.py
    order_list = []
    for filename in args.file:
        chainorder = _read_series(filename, chains_count)
        order_list.append(chainorder)
    
    #
    f = plt.figure(1)
    _plot_order(order_list,args.labels,f)
    f.savefig(args.out,format="png")
