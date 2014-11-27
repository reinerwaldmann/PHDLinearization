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


def cutPlanToTwoList(plan):
    _1=list()
    _2=list()
    for point in plan:
        _1.append(point[0])
        _2.append(point[1])
    return _1, _2



def plotPlan(plan, title=''):
    _1,_2=cutPlanToTwoList(plan)
    plt.plot(_1, _2,  'bo')
    plt.legend(loc='upper right')
    plt.title(title)
    plt.ylabel('-1')
    plt.xlabel('_2')
    plt.grid()
    plt.show()


