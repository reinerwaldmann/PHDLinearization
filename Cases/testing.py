__author__ = 'reiner'

import unittest


def compare (a,b):
    """
    returns
    1 (True), if a>b
    else 0 (False)
    Looking only at the exponent
    """
    def get_normalized_exp(x):
        """
        Получает экспоненту у нормализованного представления числа
        2 этапа: если оно имеет порядок 1 и выше, то срабатывает первая часть,
        до первого return (она уменьшает число последовательно)

        Иначе к делу подключается вторая

        """
        i=0

        while (int(x)):
            x/=10
            i+=1
        if i:
             return i
        while not (round(x,i)):
            i+=1
        return -1*i
    return get_normalized_exp(a)>get_normalized_exp(b)



#testing


class testCompareForKnownValues (unittest.TestCase):
    def testfractForward(self):
        self.assertTrue(compare(.1,.001))
    def testfractReverse(self):
        self.assertFalse(compare(.001,.1))
    def testfractEq(self):
        self.assertFalse(compare(.2,.1))

    def testforward(self):
        self.assertTrue(compare(100.2,1))

    def testreverse(self):
        self.assertFalse(compare(1, 100.2))

    def testeQ(self):
        self.assertFalse(compare(300, 100.2))

    def testmixForward(self):
        self.assertTrue(100,.0001)

    def testmixReverse(self):
        self.assertTrue(.0001, 100)



#
if __name__ == '__main__':
    unittest.main()

#     test()




