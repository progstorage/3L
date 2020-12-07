from math import *
import numpy as np
import matplotlib.pyplot as plt
import random as rnd

def Fun(x, h):
	fi1 = []
	r1 = []

	if (x == 1):
		a = 0.1
		for fi in np.arange(0, 2 * pi, h):
			fi1.append(fi)
			if ((2 * cos(2 * fi)) >= 0):
				r1.append(a * sqrt(2 * cos(2 * fi))) #r^2=2*(a^2)*cos(2*fi)
			else:
				r1.append(a * sqrt(-2 * cos(2 * fi)))

	elif (x == 2):
		for fi in np.arange(0, 2 * pi, h):
			fi1.append(fi)
			r1.append(2 - 4*sin(fi)) #r=2-4*sin(fi)

	elif (x == 3):
		a = 1.7
		for fi in np.arange(0, a * pi, h):
			fi1.append(fi)
			r1.append(sin((7 * fi) / 4)) #r=sin(7*fi/4)

	elif (x == 4):
		a = 1
		for fi in np.arange(0, 2 * pi, h):
			fi1.append(fi)
			r1.append(a * cos(3 * fi)) #r=a*cos(3*fi)

	elif (x == 5):
		for fi in np.arange(0, 2 * pi, h):
			fi1.append(fi)
			r1.append((cos(fi)**2) + 2*sin(fi)) #r=cos(fi)^2 + 2*sin(fi)

	elif (x == 6):
		y = lambda t: (2*cos(3 * t) + 3*sin(0.5 * t)) / 0.5*cos(4 * t)
		(a, b) = (-2, 2)
		_size = int((b - a) / h)
		Y = [[0] * 2 for i in range(_size)]
		i = 0
		for _x in np.arange(a, b, h):
			Y[i][0] = _x
			Y[i][1] = y(_x)
			i += 1
		return Y


	Y = [[0] * 2 for i in range(len(fi1))]

	for i in range(len(fi1)):
		Y[i][0] = (r1[i] * cos(fi1[i]))
		Y[i][1] = (r1[i] * sin(fi1[i]))
	return Y

h = 0.02 #шаг
fun_count = 6
for i in range(1, fun_count + 1):
	xy = Fun(i, h)  
	# для отрисовки полученной функции:
	plt.scatter([row[0] for row in xy], [row[1] for row in xy], 1)	
	plt.show()

	x_noise = rnd.uniform(-5, 5)
	y_noise = rnd.uniform(-5, 5)

	file_name = '../in/in_{}.txt'.format(i)
	with open(file_name, 'w+') as _in:	
		for j in range(len(xy)):
			_in.write(str(xy[j][0] + x_noise))
			_in.write(';')
			_in.write(str(xy[j][1] + y_noise))
			_in.write('\n')



















