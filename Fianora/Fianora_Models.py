__author__ = 'vasilev_is'

class AbstractModel:
    def funcf (self, x,b):
        raise NotImplementedError("funcf Please Implement this method")

    def jacf(self, x, b, y):
        raise NotImplementedError("jacf Please Implement this method")

    def hess(self, x, b, y):
        raise NotImplementedError("hess Please Implement this method")


    def cal_funcf(self):
        return self.funcf

    def cal_jacf(self):
        return self.jacf

    def cal_hess(self):
        return self.hess


class SemiconductorModel(AbstractModel):
    FT = 0.02586419


class SimpleDiodeModel(SemiconductorModel):
    pass


class AdvancedDiodeModel(SemiconductorModel):
    pass


class EMTransistorModel(SemiconductorModel):
    pass


class GPTransistorModel(SemiconductorModel):
    pass
