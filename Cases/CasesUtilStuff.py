"""
Здесь размещаются инструменты для создания и обслуживания кейсов различных

"""


__author__ = 'vasilev_is'


class IterationInfoAcceptor ():
    """
    Это класс, который ведёт себя по сути как список
    он принимает на вход любую хрень, и складирует её в список,
    тот список можно получить через переопределённые стандартные операторы
    Если задан файл, то он будет туда писать по мере прихода данных
    При начальном обявлении файл перезаписывается

    """

    def __init__(self, filename=None, verbose=False, read=False):
        self.lst_of_data=[]
        self.verbose = verbose
        self.lineCounter = 0

        if read and filename:
            with open(filename, 'r') as f:
                for line in f:
                    ll = line.split('\t')
                    try:
                        lst = eval(ll[1])
                    except:
                        lst=None

                    try:
                        dct = eval(ll[2])
                    except:
                        dct=None

                    self._smallaccept(lst, dct)
            return

        if filename:
            self.filename = filename
            with open(self.filename, 'w') as f:
                pass


    def _smallaccept (self, args, kwargs):

        if args and kwargs:
            self.lst_of_data.append((args,kwargs))
            return

        if args:
            self.lst_of_data.append(args)

        if kwargs:
            self.lst_of_data.append(kwargs)



    def accept (self, *args, **kwargs):

        if hasattr(self, 'filename') and self.filename:
                with open(self.filename, 'a') as f:
                    f.write(str(self.lineCounter)+"\t"+args.__str__()+"\t"+kwargs.__str__()+"\n")


        if self.verbose:
            print (self.lineCounter, args, kwargs)

        self.lineCounter+=1

        self._smallaccept(args, kwargs)



    def __getitem__(self, key):
        return self.lst_of_data[key]

    def __iter__(self):
        return self.lst_of_data.__iter__()
