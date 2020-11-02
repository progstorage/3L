import math
import numpy as np

def f(x):
	rez=math.atan(x)
	return rez

n=500 #количество узлов
x=[]
y=[]

h=0.1
import numpy as np
first_value=np.random.uniform(0,100) 
for i in range(n):
	x.append(round(first_value,5))
	first_value=first_value+h

for i in x:
	y.append(round(f(i),5))

for i in range(len(y)):
	print('x: ', x[i],'| y: ', y[i])


import matplotlib.pyplot as plt
plt.plot(x, y, linewidth=1.0)
plt.grid()
plt.title('(График функции)')
plt.show()


with open("array.txt", "w") as f1:
    for  line in range(len(x)):
        f1.write(str(x[line]) + ';' + str(y[line]) + '\n')
f1.close ()
