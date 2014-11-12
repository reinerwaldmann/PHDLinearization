__author__ = 'reiner'




def getderivative_outTransParamErlieFormatJAC():
    """
    Якобиан по коэфффициентам
    """
    funstr=["B2*( (1/B0)*(exp(Vbe/FT)-1)+(1/B1)*(exp(Vbc/FT)-1)) + GMIN*(Vbe/B0+Vbc/B1)", "B2*( (exp(Vbe/FT)-exp(Vbc/FT))*(1-Vbc/B3)-(1/B1)*(exp(Vbc/FT)-1))+GMIN*((Vbe-Vbc)*(1-Vbc/B3)-Vbc/B1)" ]

    resstr=""

    for i in range (0, len(funstr)):
        for ind in range (0, 4):
            #print(sympy.diff(funstr[i],'B{0}'.format(ind)))


            resstr+=sympy.diff(funstr[i],'B{0}'.format(ind)).__str__()
            resstr+="\n"
            print ('B{0}'.format(ind))

        resstr+="------------------\n"

    return resstr

