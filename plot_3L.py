import matplotlib.pyplot as plt


x = [-1]
y = [0]
i = -0.8
while i<=1:
	x.append(i)
	y.append((1-i*i)**(1/2))
	i += 0.2

l1 = []
l3 = []
with open(rf'D:\u\ConsoleApplication1\ConsoleApplication1\out.txt', 'r') as f:
	for line in f:
		l = line.strip('\n').split(',')
		l1.append([float(l[0]), float(l[1])])
		l3.append([float(l[2]), float(l[3])])

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.grid(True, which='both')
ax.set_ylim(-0.5, 2)
ax.set_xlim(-2, 2)
plt.scatter(x, y)

plt.scatter([i[0] for i in l1], [i[1] for i in l1])
plt.scatter([i[0] for i in l3], [i[1] for i in l3])

plt.gca().set_aspect('equal', adjustable='box')
plt.draw()
plt.show()