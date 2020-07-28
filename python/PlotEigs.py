import numpy as np
import matplotlib.pyplot as plt
import math

if __name__ == '__main__':
    numPlots = 10
    f = open('/tmp/eig10')
    lines = f.readlines()
    x = []
    y = []
    for line in lines:
        tokens = line.strip().split(' ')
        x.append(float(tokens[0]))
        #y.append(math.log(abs(float(tokens[1]))))
        y.append(float(tokens[1]))
    fig, ax = plt.subplots()
    ax.plot(x,y)
    #ax.grid(False,which='both')
    ax.axhline(y=0, color='k')
    ax.axvline(x=0, color='k')
    plt.show()
