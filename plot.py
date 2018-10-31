#! /home/marie/Utilities/anaconda3/bin/python

import argparse
import numpy as np
import pylab as plt
import os
import seaborn as sns


def plot_histo():
    # Number of random variables to be generated
    nsamp = 100000

    nbins = 100

    # Initialize some variables
    xav  = 0.0
    x2av = 0.0
    histo = np.zeros(nbins)
    x = np.zeros(nbins)
    xmin =-10.0
    xmax = 10.0

    # Histogram bin-width
    deltax = (xmax - xmin) / nbins


    # Generate random numbers
    for i in range (nsamp):

        # Uniform random numbers in [0,1]
        rnum1 = np.random.random_sample()
        rnum2 = np.random.random_sample()

        # Gaussian variable
        rnum3 = np.sqrt(-2.0 * np.log(rnum1)) * np.cos(2.0 * np.pi * rnum2)

        xav  = xav + rnum3
        x2av = x2av + pow(rnum3,2.0)

        #Calculate histogram here

        ibin = int((rnum3 - xmin) / deltax)
        histo[ibin] += 1.0


    # Calculate mu and sigma
    mu = xav / float(nsamp)
    sigma = np.sqrt(x2av / float(nsamp) - pow(mu,2.0))

    #Write mu and sigma
    print("mu     =", mu)
    print("sigma  =", sigma)


    # set histogram for plotting
    norm = 0.0
    for i in range(nbins):

        x[i] = xmin + (i * deltax)
        histo[i] = histo[i] / nsamp
        norm += histo[i]


    width = 0.7 * (x[1] - x[0])
    plt.bar(x, histo, color='green', align='center', width=width)
    plt.title("Gaussian Histogram")
    plt.xlabel("Value")
    plt.ylabel("Frequency")

    plt.show()


def plot(files):

    plt.figure()
    legend=[]

    for file in files:
        print('reading file')
        datas=np.loadtxt(file,  delimiter=',',unpack=True, comments='#')
        with open(file) as f:
            first_line = f.readline().replace("#", "").replace("\n", "").split(',')
        for col in range(1,len(first_line)):
            if col+1 < len(first_line) and 'err' in first_line[col+1]:
#                plt.errorbar(datas[0],datas[col], yerr=datas[col+1],ecolor='k',mfc='black',capsize=3,linewidth=1,label=first_line[col])
                plt.errorbar(datas[0],datas[col], yerr=datas[col+1],capsize=3,linewidth=1,label=first_line[col])
                legend.append(first_line[col])
            elif 'err' not in first_line[col] :
                plt.errorbar( datas[0],datas[col],label=first_line[col],linewidth=1)
                legend.append(first_line[col])
    plt.legend()


def main():

    parser = argparse.ArgumentParser(description='plot files')
    parser.add_argument('files', type=str,nargs='*',
                        help='files for ploting')
    parser.add_argument('-title', default=False,
                        help='tittle displayed on the graph')
    parser.add_argument('-xaxis', default=False,
                        help='xaxis,yaxis')
    parser.add_argument('-yaxis', default=False,
                            help='xaxis,yaxis')
    args = parser.parse_args()


    flatui = [ "#3498db", '#ff3333', "#59b69c",  "#34495e","#b69c59", "#2ecc71" ,"#9b59b6" , "#95a5a6"]
    sns.set_palette(flatui)

    plot(args.files)
### axis and stuffs

    plt.title(args.title)
    plt.xlabel(args.xaxis)
    plt.ylabel(args.yaxis)



    plt.savefig('graph.png', bbox_inches='tight')

if __name__ == "__main__":
    main()
