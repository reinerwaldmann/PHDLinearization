__author__ = 'SW'

import matplotlib.pyplot as plt
import numpy as np

def plotPlanAndMeas2D(measdata, title=''):

    planplot1=[x['x'][0] for x in measdata]
    measplot1=[x['y'][0] for x in measdata]
    plt.plot(planplot1, measplot1,  'ro')
    plt.title(title)
    plt.ylabel('value')
    plt.xlabel('point')
    plt.grid()
    plt.show()





def plotSkGraph(gknu, title=''):
    #plotting Sk graph
    rng=np.arange(0,len(gknu[3]))
    plt.plot(rng , gknu[3], label='Sk drop')
    plt.legend(loc='upper right')
    plt.title(title)
    plt.ylabel('Sk')
    plt.xlabel('Interation')
    plt.grid()
    plt.show()
