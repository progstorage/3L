import matplotlib.pyplot as plt
import numpy as np

def data():
	'''x = [-1]
	y = [0]
	i = -0.8
	while i<=1:
		x.append(i)
		y.append((1-i*i)**(1/2))
		i += 0.2
	with open('data.txt', 'w') as f:
		for i in range(len(x)):
			f.write(str(x[i])+';' + str(y[i])+'\n')'''
	x = []
	y = []
	with open('data.txt', 'r') as f:
		for line in f:
			x.append(float(line.split(';')[0]))
			y.append(float(line.split(';')[1]))
	return (x, y)


def plot_g(x, y):
	l1 = []
	l3 = []
	with open('out.txt', 'r') as f:
		for line in f:
			l = line.strip('\n').split(',')
			l1.append([float(l[0]), float(l[1])])
			l3.append([float(l[2]), float(l[3])])

	a = []
	with open('coef.txt', 'r') as f:
		for line in f:
			a.append(float(line))

	fig = plt.figure()
	ax = fig.add_subplot(1, 1, 1)
	ax.grid(True, which='both')
	ax.set_ylim(min(y)-1, max(y)+1)
	ax.set_xlim(min(x)-1, max(x)+1)
	plt.scatter(x, y)
	plt.scatter([i[0] for i in l1], [i[1] for i in l1])
	plt.scatter([i[0] for i in l3], [i[1] for i in l3])

	plt.plot(x,y)
	plt.plot([i[0] for i in l1], [i[1] for i in l1])
	plt.plot([i[0] for i in l3], [i[1] for i in l3])

	f = lambda x, y: a[0]+a[1]*x+a[2]*y+a[3]*x*x+a[4]*x*y+a[5]*y*y
	delta = 0.025
	x_range = np.arange(min(y)-1, max(y)+1, delta)
	y_range = np.arange(min(x)-1, max(x)+1, delta)
	X, Y = np.meshgrid(x_range, y_range)
	F = f(X, Y)
	plt.contour(X, Y, F, [0], colors = ['red'])

	plt.gca().set_aspect('equal', adjustable='box')
	plt.draw()
	plt.show()
	

if __name__ == '__main__':
	x, y = data();
	plot_g(x, y);