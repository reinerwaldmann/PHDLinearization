__author__ = 'vasilev_is'

import pickle

import matplotlib.pyplot as plt
from Cases.CasesUtilStuff import IterationInfoAcceptor
import numpy as np
import math

def norm1 (b1,b2):
    absdif = [math.fabs(b11-b22) for b11,b22 in zip(b1,b2)]
    return max(absdif)

def norm2 (b1,b2):
    return np.linalg.norm(np.array(b1)-np.array(b2))

def norm3 (b1,b2):
    absdif = [math.fabs(b11-b22) for b11,b22 in zip(b1,b2)]
    return sum(absdif)



datafile = 'resfiles/resdump.dat'

with open(datafile, 'rb') as f:
    rrlist = pickle.load(f)

# for case in rrlist:
#     print (case.ll())


rrlist_filtered = [case for case in rrlist if case.ll()<20]

print ('dif', len(rrlist)-len(rrlist_filtered))



#для каждой итерации создаём гистограмму рассеяния изменений

for it in range (10):
    #it=1
    dflist=[]  # список изменений на первой итерации по кейсам
    for case in rrlist_filtered:
        #dflist.append(norm3(case[0]['b'], case[1]['b']))
        try:
            dflist.append(case[it]['b'][2]- case[it+1]['b'][2])
        except:
            pass

    fig, ax = plt.subplots( nrows=1, ncols=1 )  # create figure & 1 axis
    ax.hist(dflist,30)

    #это у нас графики значений компонента вектора по итерациям.
    fig.savefig ('resfiles/changes_hists_b2/img_{0}.png'.format(it))
    plt.close(fig)


    #bonelist = [zp['b'][2] for zp in case]
    #plt.plot(bonelist)
#plt.hist(list(filter(lambda x: x<1.1, dflist)), 50)
#plt.show()



