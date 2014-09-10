#!/usr/bin/python


def find_nearest(vector, value):
    '''
    Find value in \'vector\' nearest to \'value\'
    '''
    import numpy as np
    return np.abs(vector - value).argmin()


def avg1(a):
    return (a[:-1]+a[1:])/2
