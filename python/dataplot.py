#!/usr/bin/env python
from astropy.io import ascii
import matplotlib.pyplot as plt

data=ascii.read('../data/RXJ1131_Tewes2013.rdb',data_start=2)
x=data['mhjd']
y=data['mag_A']
plt.scatter(x,y,c='r')
y=data['mag_B']
plt.scatter(x,y,c='b')
plt.show()


