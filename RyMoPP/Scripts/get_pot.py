#!/usr/bin/python
from decimal import *
from numpy import *
import glob as gb
getcontext().prec = 6
files = gb.glob('EigenValue*.txt')
start = 22.00000
names = arange(start,40.4000000,0.2)
data = []
for i in names:
        eigen = loadtxt('EigenValue'+ str('{0:.6f}'.format(i)) + '.txt')
        if isinstance(eigen,float):
                eigenabs = [abs(eigen)]
                eigen = [eigen]
        else:
                eigenabs = abs(eigen).tolist()
                eigen = eigen.tolist()
        eigen = [i] + eigen
        data.append(eigen)
savetxt('data.txt', array(data))
