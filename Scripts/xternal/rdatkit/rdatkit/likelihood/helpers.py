from numpy import *


def normalize(values):
    a = array(values)
    m = a[logical_not(isnan(a))].mean()
    a = a*(1./m)
    a[a < 0.01] = NaN
    return a


