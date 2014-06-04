__author__ = 'vasilev_is'


import numpy as np

def generate (funcstrdict, xvectorlistsdict, spreadvarslist=None, Vx=None, nvolx=1, yvectordispsdict=None, nvoly=1, outfilename=""):
    """
Функция способна в широких пределах моделировать экспериментальные данные, получаемые в процессе исследования
сущности, задаваемой многооткликовой моделью.
Описание прилагается

Входные параметры:
словарь входных функций  - выходная переменная к строке функции
словарь значений входных переменных  - входная переменная к списку значений (возможно, задаваемому как arange, range или каким-либо ещё образом). Если вх. перем. имеют разброс, то это список матожиданий
список переменных с разбросом, в строковом представлении
V, ковариационная матрица переменных с разбросом. Диагональная V говорит о наличии дисперсии, но отсутствии ковариации
nvolx, размер выборки входных переменных, получаемых с разбросом, для каждого матожидания
словарь дисперсий выходных переменных - выходная переменная - дисперсия.
nvoly, размер выборки выходных переменых, получаемых с разбросом, для каждого матожидания
имя выходного файла, если имеется потребность в записи входных и выходных переменных в выходной файл


Выходные значения:
список словарей выходная переменная - значение
опционально файл, где записаны входные переменные - выходные переменные. Столбцы озаглавливаются.

    """
    res=dict()


    evalstr=""
    counter=0
    for key in xvectorlistsdict.keys():
        for i in range (0, counter):
            evalstr+="\t"
        evalstr+="for "+key+" in "+xvectorlistsdict[key].__str__()+":\n"
        counter+=1







    resdict=dict()

    for key in funcstrdict.keys():

        for i in range (0, counter):
            evalstr+="\t"
        resdict[key]=list()
        evalstr+="resdict[\""+key+"\"].append("+funcstrdict[key]+")\n"


    exec(evalstr, {"u1":0, "r3":0, "r2":0, "r1":0 }, locals() )

    print (resdict)


    print (evalstr)






    #for x in xvectorlistsdict["r1"]:
        #print (x)



    return res








funcstrdict= {"y1":"u1* (r2+r3)/(r1+r2+r3)", "y2":"u1* r3/(r1+r2+r3)"}
xvectorlistsdict = {"u1":[100],  "r1":list(np.arange(0,0.5,0.1)), "r2":list(range(0,3,1)), "r3":list(range (20,30,5))}



generate (funcstrdict, xvectorlistsdict, spreadvarslist=None, Vx=None, nvolx=1, yvectordispsdict=None, nvoly=1, outfilename="")







#
# def test (inp):
#     for j in list (inp):
#         print (j)


#test (np.arange(0,1,0.1))



