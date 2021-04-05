import random as rand
import matplotlib.pyplot as plt
import numpy as np

for i in range(20, 100001):
    a0 = rand.uniform(-10000, 10000)
    a1 = rand.uniform(-10000, 10000)
    a2 = rand.uniform(-10000, 10000)
    a3 = rand.uniform(-10000, 10000)
    a4 = rand.uniform(-10000, 10000)
    a5 = rand.uniform(-10000, 10000)
    #a6 = rand.uniform(-10000, 10000)
    #a7 = rand.uniform(-10000, 10000)
    
    b0 = rand.uniform(-10000, 10000)
    b1 = rand.uniform(-10000, 10000)
    b2 = rand.uniform(-10000, 10000)
    b3 = rand.uniform(-10000, 10000)
    b4 = rand.uniform(-10000, 10000)
    b5 = rand.uniform(-10000, 10000)
    #b6 = rand.uniform(-10000, 10000)
    #b7 = rand.uniform(-10000, 10000)
    
    A = lambda x: a0 + a1*x + a2*x**2 + a3*x**3 + a4*x**4 + a5*x**5 #+ a6*x**6 + a7*x**7
    B = lambda x: b0 + b1*x + b2*x**2 + b3*x**3 + b4*x**4 + b5*x**5 #+ b6*x**6 + b7*x**7
    
    arr = []
    for x in np.arange(-150, 150, 1):   # 300 points for each test
        arr.append( [x, A(x) / B(x)] )
    x_noise = rand.uniform(-100, 100) / 1000
    y_noise = rand.uniform(-100, 100) / 1000

    file_name = '../in/in_{}.txt'.format(i)
    with open(file_name, 'w+') as _in:	
        for j in range(len(arr)):
            _in.write(str(arr[j][0] + x_noise))
            _in.write(';')
            _in.write(str(arr[j][1] + y_noise))
            _in.write('\n')
