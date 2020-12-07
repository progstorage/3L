# отрисовывает все массивы из /tests/in

import os
import matplotlib.pyplot as plt

def draw(file):
	x = []
	y = []
	with open(file, 'r+') as _in:
		for line in _in:
			x.append(float(line.split(';')[0]))
			y.append(float(line.split(';')[1]))
		plt.scatter(x, y)
		plt.plot(x, y)
		plt.show()


dir_name = 'in/'
for _, _, files in os.walk(dir_name):
	i = 0
	for file in files:
		i += 1
		plt.figure(i)
		draw(dir_name + file)
