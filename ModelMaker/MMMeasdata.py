__author__ = 'vasilev_is'

import csv


#TODO привести всё в соответствии с PEP8 и PEP257
class MMMeasdata:
    """
    Содержит в себе данные из датасета, принимаемые на вход методом init (ему даётся имя файла)
    init предполагает, что все данные-числовые, притом float
    Далее, методом setAccessList запиываем идентификаторы переменных, которые станут членами вектора x и вектора y
    После этого объект из себя представляет measdata - можно делать доступ по квадратным скобкам, можно итеративным процессом
    """


    def __init__(self, filename):
        """
        этот метод должен считывать данные откуда-то, лучше всего с csv файла
        """
        with open(filename, 'rt') as f:
            reader=csv.reader(f)

            self.data=list()
            rnum=0
            for r in reader:
                if rnum==0:
                    self.readableNames = r
                elif rnum==1:
                    self.ids = r
                else:
                    self.data.append([float (rr) for rr in r] )
                rnum+=1

        self.inlist=[]
        self.outlist=[]


    def __len__(self):
        return len(self.data)

    def setAccessList(self, inlist, outlist):
        """
        Устанавливаем, кто будет членами  вектора x, кто будет членами вектора y
        :param inlist:
        :param outlist:
        :return:
        """
        self.inlist = inlist
        self.outlist = outlist

    def __getitem__(self, item):
        """
        Эмулирует из себя measdata, то есть список словарей с ключами x и y
        :param item:
        :return:
        """
        #проверка на установленность атрибутов inlist, outlist

        if self.inlist and self.outlist:
            invarsindexes = [self.ids.index(x) for x in self.inlist]  #индексы требуемых входных переменных
            outvarsindexes = [self.ids.index(x) for x in self.outlist]  #индексы требуемых выходных переменных

            xvector = [self.data[item][i] for i in invarsindexes]
            yvector = [self.data[item][i] for i in outvarsindexes]

            return {'x':xvector, 'y':yvector}

        else:
            raise AttributeError('Inlist or/and Outlist not set, thus method cannot work')

    def showReadableNames(self,shouldBeFromModel=False, shouldBeIds=False):
        """
        Показывает читаемые названия полей
        :param shouldBeFromModel: если тру, тогда только те, которые в модели (in, out)
        :return:
        """
        if shouldBeIds: pl=self.ids
        else: pl=self.readableNames

        print("Variabnles in model description")
        if shouldBeFromModel:
            print("Input Variables: ")
            print([pl[self.ids.index(i)] for i in self.inlist ])
            print("Output Variables: ")
            print([pl[self.ids.index(i)] for i in self.outlist ])
        else:
            print (pl)


    def getCut(self,id_):
        """
        получает срез  (список переменной) по идентификатору
        :param id:
        :return:
        """
        idx = self.ids.index(id_)
        return [lst[idx] for lst in self.data]

    def getXarray(self):
        """
        :return:
        """
        return [v['x'] for v in self]
    def getY(self):
        """
        :return:
        """
        return [v['y'][0] for v in self]





def test():
    measdata = MMMeasdata('test.csv')

    measdata.inlist=['in1', 'in2']
    measdata.outlist=['out1']

    for m in measdata:
        print (m)


#test()