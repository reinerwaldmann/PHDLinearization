__author__ = 'reiner'

import unittest
from  Cases.MantissaEndDifferences import compare_exp




#testing
#
#
class testCompareForKnownValues (unittest.TestCase):
    def testEqFract(self):
        self.assertTrue(compare_exp(.00001, .00008))

    def testEqFract2(self):
        self.assertTrue(compare_exp(.1, .1))

    def testEq(self):
        self.assertTrue(compare_exp(100,200))

    def testEqThr(self):
        self.assertTrue(compare_exp(.1,.01,2))

    def testEqFlsThr(self):
        self.assertFalse(compare_exp(.1,.001,2))

    def testEqFlsThr1(self):
        self.assertFalse(compare_exp(.001,.1,2))

    def testEqFls(self):
        self.assertFalse(compare_exp(.001,.1))



    # def testfractForward(self):
    #     self.assertTrue(compare(.1,.001))
    # def testfractReverse(self):
    #     self.assertFalse(compare(.001,.1))
    # def testfractEq(self):
    #     self.assertFalse(compare(.2,.1))
    #
    # def testforward(self):
    #     self.assertTrue(compare(100.2,1))
    #
    # def testreverse(self):
    #     self.assertFalse(compare(1, 100.2))
    #
    # def testeQ(self):
    #     self.assertFalse(compare(300, 100.2))
    #
    # def testmixForward(self):
    #     self.assertTrue(100,.0001)
    #
    # def testmixReverse(self):
    #     self.assertTrue(.0001, 100)



#
if __name__ == '__main__':
    unittest.main()
#     compare(.1,.5)
#
#     test()




