from math import *
import numpy as np
import matplotlib.pyplot as plt


y = lambda x: (2*cos(3*x)+3*sin(0.5*x)) / 0.5*cos(4*x)

x1 = []
y1 = []
for x in np.arange(-2, 2, 0.025):	# 160 points
	x1.append(x)
	y1.append(y(x))
	
with open('../in/in_1.txt', 'w+') as _in:	
	for i in range(len(x1)):
		_in.write(str(x1[i]))
		_in.write(';')
		_in.write(str(y1[i]))
		_in.write('\n')
		
# для отрисовки полученной функции:
#plt.plot(x1, y1)	
#plt.show()