__author__ = 'vasilev_is'
import os
#функция, которая считывает имя файла в measdata

def readMeasdata (infilename):
    rs=list()

    with open (infilename) as f:
        lines = f.readlines()

    for line in lines:
        try:
            bl=line.split(',')
            x=[float(bl[0].strip())]
            y=[float(bl[1].strip())]
            rs.append( {'x':x, 'y':y } )


        except:
            pass
    if not rs:
        print ('empty'+infilename)

    rs.sort(key=lambda x: x['x']) #сортировочка для тех случаев, когда план не упорядочен


    return rs


def grfiles (folder):
    rs=list() #то будет список списков
    for file in os.listdir(folder):
        if file.endswith('.txt'):
            rs.append((readMeasdata(folder+'/'+file), file) )
    return rs


def test ():
    #print (readMeasdata('D:\PHd\RealDiode\\1N4004\Reports\\adv\\01.txt'))
    print (grfiles('D:\PHd\RealDiode\\1N4004\Reports\\adv\\'))


