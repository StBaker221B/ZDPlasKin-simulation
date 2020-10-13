
import numpy as np 
from matplotlib import pyplot as plt 
import matplotlib 


def readdata(datafile, indexfile = 0):
    
    qtfile = False
    if(indexfile):
        index = readlist(indexfile)
        qtfile = True
    
    f = open(datafile, 'r')
    data = []

    if(qtfile):
        title = f.readline()
        title = title.split()
        data.append(index)
        for i in title:
            data.append([])
    else:
        line = f.readline().split()
        index = {}
        data.append(index)
        for i in range(len(line)):
            index[i] = line[i]
            data.append([])
    # title = f.readline()
    # title = title.split()
    # print(len(title))
    # print(title)

    # for i in title:
    #     data.append([])
    line = f.readline()
    while(line):
        line = line.split()
        for i in range(len(line)):
            data[i+1].append(float(line[i]))
        line = f.readline()
    f.close()

    return data

def readlist(filename):
    f = open(filename, 'r')
    index = {}
    line = f.readline()
    while(line):
        line = line.split()
        index[int(line[0])] = line[1:]
        line = f.readline()
    f.close()
    # print(index)
    return index
    
def plotdata(data, index=[1]):
    time = data[1]
    var = []
    plt.figure(figsize=(16,9)) 
    item = ''
    for i in index:
        # var.append(data[i+1])
        print(i)
        var = data[i+1]
        item = ''
        for s in data[0][i]:
            item = item + '_' + str(s)
        # print(type(item))
        plt.plot(time, var, label = data[0][i] )
    
    fs = 20
    plt.xlim(xmin=0)
    # plt.xscale('log')
    plt.xlabel('Time [s]', fontsize = fs)
    # plt.rcParams['xtick.direction'] = 'in'
    plt.ylim( ymin=1.0e0 )
    # plt.ylim( ymin=1.0e17 )
    # plt.yscale('log')
    plt.ylabel('Number density [cm^-3]', fontsize = fs)
    # plt.rcParams['ytick.direction'] = 'in'
    plt.tick_params(top='on',bottom='on',left='on',right='on', labelsize = fs)
    plt.tick_params(direction='in')

    # plt.margins(0,0)
    # plt.axis('off')

    num1 = 1.05
    num2 = 0
    num3 = 3
    num4 = 0

    plt.legend(bbox_to_anchor=(num1, num2), loc=num3, borderaxespad=num4, fontsize = fs/2)
    # plt.legend( fontsize = fs)

    # plt.title(data[0][i], fontsize = 2*fs)
    plt.title('Time evolution', fontsize = 2*fs)

    # plt.subplot(122)
    # plt.plot(x, y3) 
    # plt.xscale('log')
    # plt.yscale('log')
    # plt.title('Te')
    # plt.plot(y, y,'ro')
    # plt.ylabel()
    # plt.xlabel('Time')
    # plt.legend()

    # figname = 'plot_' + data[0][i] + '.jpg'
    # figname = 'plot' + item + '.png'
    # figname = 'plot reactions.png'

    figname = 'plot_conditions.png'

    plt.savefig(figname, bbox_inches='tight')

def plotslice(data, time, index = [1]):
    y = []
    x = range(1, len(index)+1)
    for t in range(len(data[1])):
        # print(data[1][t])
        # print(type(data[1][t]))
        if( data[1][t] >= time ):
            for i in index:
                y.append(data[i+1][t])
            break
    # print(x)    
    # print(y)
    plt.figure(figsize=(16,9))
    plt.plot(x, y, 'ro-')

    fs = 20
    # plt.xlim( xmin = 1e-12, xmax = 2e-5)
    # plt.xscale('log')
    plt.xlabel('CO2 vibrational state', fontsize = fs)
    plt.xticks(x)
    plt.ylim( ymin = 1.0e8, ymax = 1e20 )
    plt.yscale('log')
    plt.ylabel('Number density [cm^-3]', fontsize = fs)
    plt.tick_params(top='on',bottom='on',left='on',right='on', labelsize = fs)
    plt.tick_params(direction='in')

    num1 = 1.05
    num2 = 0
    num3 = 3
    num4 = 0

    # plt.legend(bbox_to_anchor=(num1, num2), loc=num3, borderaxespad=num4, fontsize = fs/2)

    plt.title('Time slice', fontsize = 2*fs)

    figname = 'plot_' + str(time) + 's.png'
    plt.savefig(figname, bbox_inches='tight')




if __name__ == '__main__':

    # indexfile = 'qt_species_list.txt'
    # datafile = 'qt_densities.txt'
    
    indexfile = 'qt_conditions_list.txt'
    datafile = 'qt_conditions.txt'

    # indexfile = 'qt_reactions_list.txt'
    # datafile = 'qt_rates.txt'

    data = readdata(datafile, indexfile)
    # plotdata(data,range( 6, 27 ))
    # plotslice(data, 1e-9, range(6, 27))
    # plotslice(data, 1e-10, range(6, 27))
    # plotslice(data, 1e-8, range(6, 27))
    # plotslice(data, 1e-8+10e-9, range(6, 27))
    # plotslice(data, 1e-8+30e-9, range(6, 27))
    # plotslice(data, 1e-8+100e-9, range(6, 27))
    # plotslice(data, 1e-8+1e-6, range(6, 27))
    # plotslice(data, 1e-8+6e-9, range(6, 27))

    # plotdata(data,range( 668, 688 ))
    # plotdata(data, range( 617, 642 ))

    # plotdata(data, [1, 30, 60])
    # plotdata(data, [ 1, 30, 59, 60, 67, 70, 71])
    # plotdata(data, [ 1, 2, 3, 4, 5 ])
    # plotdata(data, [ 6, 7, 8, 9 ])

    plotdata(data, [1, 4, 5])
