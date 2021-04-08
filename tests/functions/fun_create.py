from math import *
import cmath
#from mpmath import *
import numpy as np
import matplotlib.pyplot as plt
import random as rnd

def Fun(x, h):
	fi1 = []
	r1 = []

	x_ = []
	y_ = []

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

	elif (x == 7):
		for i in np.arange(-1.8,1.8,h):
			x_.append(i)
			y_.append(-(1/2)+sin(i**2))

			if (i>-3/2) and (i<3/2):
				x_.append(i)
				y_.append(i**4 - 2*(i**2))

			if (i>-1/2) and (i<1/2):
				x_.append(i)
				y_.append(1/5+(i**2)*10)

		Y = [[0] * 2 for i in range(len(x_))]
		for i in range(len(x_)):
			Y[i][0] = x_[i]
			Y[i][1] = y_[i]
		return Y

	if (x == 8):
		for i in np.arange(0,5,h):
			if (i!=3):
				x_.append(i)
				y_.append((2/(i-3)+2))
				x_.append(-i)
				y_.append((2/(i-3)+2))
		Y = [[0] * 2 for i in range(len(x_))]
		for i in range(len(x_)):
			Y[i][0] = x_[i]
			Y[i][1] = y_[i]
		return Y

	if (x == 9):
		for i in np.arange(-1,1,h):
			x_.append(i)
			y_.append(-1+(i**5))

			x_.append(i)
			y_.append(10*tan(i**2))

			x_.append(i)
			y_.append(10*tan(i**3 + 1/5))
		Y = [[0] * 2 for i in range(len(x_))]
		for i in range(len(x_)):
			Y[i][0] = x_[i]
			Y[i][1] = y_[i]
		return Y

	if (x == 10):
		for i in np.arange(-1,1,h):
			x_.append(i)
			y_.append(1/2-(2*i**2))

			x_.append(i)
			y_.append(i)

			x_.append(i)
			y_.append(-i)
		Y = [[0] * 2 for i in range(len(x_))]
		for i in range(len(x_)):
			Y[i][0] = x_[i]
			Y[i][1] = y_[i]
		return Y

	if (x == 11):
		for i in np.arange(-4.6,4.6,0.1):
			x_.append(i)
			y_.append(sec(1.5*i))

			x_.append(i)
			y_.append(1.5*cos(i))

		Y = [[0] * 2 for i in range(len(x_))]
		for i in range(len(x_)):
			Y[i][0] = x_[i]
			Y[i][1] = y_[i]
		return Y

	if (x == 12):
		for i in np.arange(-2,2,h):
			x_.append(i)
			y_.append(sqrt(4-(i**2)))

			x_.append(i)
			y_.append(-sqrt(4-(i**2)))

			x_.append(i)
			y_.append((i**2)-2)

			x_.append(i)
			y_.append(-(i**2)+2)
		Y = [[0] * 2 for i in range(len(x_))]
		for i in range(len(x_)):
			Y[i][0] = x_[i]
			Y[i][1] = y_[i]
		return Y

	if (x == 13):
		for i in np.arange(-0.5,1.5,h):
			x_.append(i)
			y_.append((i**5)+0.5*i-3*(i**2))

			x_.append(i)
			y_.append(1.5*(i**2)-2*i)

		Y = [[0] * 2 for i in range(len(x_))]
		for i in range(len(x_)):
			Y[i][0] = x_[i]
			Y[i][1] = y_[i]
		return Y

	if (x == 14):
		for i in np.arange(-pi/3,pi/3,h):

			x_.append(i)
			y_.append(0.5*sin(3*i))

			x_.append(i)
			y_.append(-0.5*sin(3*i))

			x_.append(i)
			y_.append(sin(i)**2)

		Y = [[0] * 2 for i in range(len(x_))]
		for i in range(len(x_)):
			Y[i][0] = x_[i]
			Y[i][1] = y_[i]
		return Y

	if (x == 15):
		for i in np.arange(-0.5,0.5,h):

			if (i<=0):
				x_.append(i)
				y_.append(2*i+1)

			if (i>=0):
				x_.append(i)
				y_.append(-2*i+1)

			x_.append(i)
			y_.append(sqrt(i+0.5))

			x_.append(i)
			y_.append(sqrt(-i+0.5))

		Y = [[0] * 2 for i in range(len(x_))]
		for i in range(len(x_)):
			Y[i][0] = x_[i]
			Y[i][1] = y_[i]
		return Y

	Y = [[0] * 2 for i in range(len(fi1))]

	for i in range(len(fi1)):
		Y[i][0] = (r1[i] * cos(fi1[i]))
		Y[i][1] = (r1[i] * sin(fi1[i]))

	return Y


def generate_polinom(power, _x):
    coef = [rnd.uniform(-1000, 1000)]
    for i in range(1, power):
        coef.append( rnd.uniform(-1000, 1000)*(_x**i) + coef[i - 1] )
    return coef
    

def generate_test(power, y):
    '''
    a_[n](x)*y^[n](x) + a_[n-1](x)*y^[n-1](x) + ... + a_[0](x) = 0,
    where:
        a_[i] - polynomial functions of x;
        y(x) - lambda
    '''
    _y = []
    _x = np.arange(-10, 10, 0.1)    # 200 points
    for t in _x:
        a = generate_polinom(power, t)
        sum = a[0]
        for i in range(1, power):
            y_noise = rnd.uniform(-100, 100) / 1000
            sum += a[i] * y(t) + y_noise
        _y.append(sum)
    
    return [_x, _y]


def run_2(x):
    if x == 0:
        power = 3
        y = lambda x: sin(x)**2
        [arr_x, arr_y] = generate_test(power, y)
         
    if x == 1:
        power = 4
        y = lambda x: cos(x)**3
        [arr_x, arr_y] = generate_test(power, y)
         
    if x == 2:
        power = 5
        y = lambda x: exp(x)**2
        [arr_x, arr_y] = generate_test(power, y)
         
    if x == 3:
        power = 6
        y = lambda x: sin(x)
        [arr_x, arr_y] = generate_test(power, y)
         
    if x == 4:
        power = 3
        y = lambda x: cos(x)**3
        [arr_x, arr_y] = generate_test(power, y)
         
    if x == 5:
        power = 2
        y = lambda x: cmath.sqrt(x)
        [arr_x, arr_y] = generate_test(power, y)
         
    if x == 6:
        power = 4
        y = lambda x: cmath.sqrt(x) * cos(x)**3
        [arr_x, arr_y] = generate_test(power, y)
             
    if x == 7:
        power = 5
        y = lambda x: cmath.sqrt(x) * sin(x)**3
        [arr_x, arr_y] = generate_test(power, y)
         
    if x == 8:
        power = 3
        y = lambda x: cmath.sqrt(x) * tan(x)**3
        [arr_x, arr_y] = generate_test(power, y)
         
    if x == 9:
        power = 2
        y = lambda x: cosh(x)**4
        [arr_x, arr_y] = generate_test(power, y)
         

    if x == 10:
        power = 6
        y = lambda x: cosh(x)**2
        [arr_x, arr_y] = generate_test(power, y)
        
    if x == 11:
        power = 5
        y = lambda x: cosh(x)**2 * tan(x)**3
        [arr_x, arr_y] = generate_test(power, y)
        
    if x == 12:
        power = 3
        y = lambda x: asinh(x)**4
        [arr_x, arr_y] = generate_test(power, y)
        
    if x == 13:
        power = 2
        y = lambda x: sin(x)**(1/2)
        [arr_x, arr_y] = generate_test(power, y)
        
    if x == 14:
        power = 5
        y = lambda x: sin(x) * cos(x)**2
        [arr_x, arr_y] = generate_test(power, y)
            
    if x == 15:
        power = 4
        y = lambda x: sin(x)**2 * cos(x)**3
        [arr_x, arr_y] = generate_test(power, y)
            
    if x == 16:
        power = 3
        y = lambda x: sin(x)**4 / tan(x)
        [arr_x, arr_y] = generate_test(power, y)
        
    if x == 17:
        power = 2
        y = lambda x: cos(x)**2 * sin(x)**4 / tan(x)
        [arr_x, arr_y] = generate_test(power, y)
        
    if x == 18:
        power = 3
        y = lambda x: sin(x) * atan(x) + cos(x)**(1/2)
        [arr_x, arr_y] = generate_test(power, y)
        
    if x == 19:
        power = 6
        y = lambda x: sin(x)**3 *cmath.sqrt(x)
        [arr_x, arr_y] = generate_test(power, y)
        

    if x == 20:
        power = 5
        y = lambda x: sin(x) * cos(x)**2 * cmath.sqrt(x)
        [arr_x, arr_y] = generate_test(power, y)
        
    if x == 21:
        power = 3
        y = lambda x: sin(x)**3 * cos(x)**2
        [arr_x, arr_y] = generate_test(power, y)
        
    if x == 22:
        power = 4
        y = lambda x: sin(x) * cosh(x)**3
        [arr_x, arr_y] = generate_test(power, y)
            
    if x == 23:
        power = 4
        y = lambda x: sin(x)**4 * exp(x)**2
        [arr_x, arr_y] = generate_test(power, y)
        
    if x == 24:
        power = 6
        y = lambda x: sin(x)**(1/3) * cmath.sqrt(x)
        [arr_x, arr_y] = generate_test(power, y)
        
    if x == 25:
        power = 7
        y = lambda x: sin(x)**3 * cmath.sqrt(x)**3
        [arr_x, arr_y] = generate_test(power, y)
        
    if x == 26:
        power = 6
        y = lambda x: sin(x) * cmath.sqrt(x)
        [arr_x, arr_y] = generate_test(power, y)
        
    if x == 27:
        power = 3
        y = lambda x: exp(x) / cmath.sqrt(x)
        [arr_x, arr_y] = generate_test(power, y)
        
    if x == 28:
        power = 4
        y = lambda x: exp(x)
        [arr_x, arr_y] = generate_test(power, y)
        
    if x == 29:
        power = 3
        y = lambda x: sin(x)
        [arr_x, arr_y] = generate_test(power, y)
        

    if x == 30:
        power = 6
        y = lambda x: exp(x) * sin(x)**4 / tan(x)
        [arr_x, arr_y] = generate_test(power, y)
            
    if x == 31:
        power = 7
        y = lambda x: exp(x) / sin(x) * cos(x)
        [arr_x, arr_y] = generate_test(power, y)
         
    if x == 32:
        power = 4
        y = lambda x: exp(x)**2 * atanh(x)
        [arr_x, arr_y] = generate_test(power, y)
        
    if x == 33:
        power = 7
        y = lambda x: exp(x)**(1/2) /cmath.sqrt(x)
        [arr_x, arr_y] = generate_test(power, y)
        
    if x == 34:
        power = 5
        y = lambda x: exp(x) * atan(x)**3
        [arr_x, arr_y] = generate_test(power, y)
        
    if x == 35:
        power = 4
        y = lambda x: exp(x) * cosh(x)
        [arr_x, arr_y] = generate_test(power, y)
        
    if x == 36:
        power = 6
        y = lambda x: exp(x) * cmath.sqrt(x)**2
        [arr_x, arr_y] = generate_test(power, y)
        
    if x == 37:
        power = 2
        y = lambda x: cos(x) * sin(x) * cmath.sqrt(x)**3
        [arr_x, arr_y] = generate_test(power, y)
        
    if x == 38:
        power = 5
        y = lambda x: cos(x) * cmath.sqrt(x)**3
        [arr_x, arr_y] = generate_test(power, y)
            
    if x == 39:
        power = 3
        y = lambda x: cos(x) * asinh(x)
        [arr_x, arr_y] = generate_test(power, y)
        

    if x == 40:
        power = 6
        y = lambda x:cmath.sqrt(x) * tan(x)**3 * tan(x)**2
        [arr_x, arr_y] = generate_test(power, y)
        
    if x == 41:
        power = 7
        y = lambda x:cmath.sqrt(x) * sin(x)**2 * tan(x)
        [arr_x, arr_y] = generate_test(power, y)
        
    if x == 42:
        power = 7
        y = lambda x:cmath.sqrt(x) * sin(x) * atan(x) + cos(x)**(1/2)
        [arr_x, arr_y] = generate_test(power, y)
        
    if x == 43:
        power = 5
        y = lambda x:cmath.sqrt(x) * sin(x)**3 + cos(x)**(1/2)
        [arr_x, arr_y] = generate_test(power, y)
        
    if x == 44:
        power = 5
        y = lambda x:cmath.sqrt(x) * cos(x)**2 * tan(x)**2
        [arr_x, arr_y] = generate_test(power, y)
        
    if x == 45:
        power = 4
        y = lambda x: tan(x) * sin(x)**3 - cos(x)**2
        [arr_x, arr_y] = generate_test(power, y)
        
    if x == 46:
        power = 3
        y = lambda x: tan(x) * sin(x) / cos(x)**2 - tan(x)
        [arr_x, arr_y] = generate_test(power, y)
            
    if x == 47:
        power = 5
        y = lambda x: tan(x)
        [arr_x, arr_y] = generate_test(power, y)
         
    if x == 48:
        power = 4
        y = lambda x: tan(x) * sin(x) / cos(x)**2 * tan(x)
        [arr_x, arr_y] = generate_test(power, y)
        
    if x == 49:
        power = 6
        y = lambda x: tan(x)**(1/2) * sin(x) * exp(x)**2
        [arr_x, arr_y] = generate_test(power, y)
            
    if x == 50:
        power = 7
        y = lambda x: tan(x) / sin(x) * cmath.sqrt(x)**3
        [arr_x, arr_y] = generate_test(power, y)
      
        
    plt.plot(arr_x, arr_y)
    plt.show()
    
    return [arr_x, arr_y]
    
    

def main():
    '''h = 0.02 #шаг
    fun_count = 15
    for i in range(10, fun_count + 1):
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
    '''
    # new tests
    fun_count_2 = 51
    for i in range(1, fun_count_2):
        [arr_x, arr_y] = run_2(i)
        if i < 10:
            file_name = '../in_2/in_00{}.txt'.format(i)
        elif i >= 10 and i < 100:
            file_name = '../in_2/in_0{}.txt'.format(i)
        else:   # i > 100
            file_name = '../in_2/in_{}.txt'.format(i)
        
        with open(file_name, 'w+') as _in:
            for j in range(len(arr_x)):
                _in.write(str(arr_x[j]))
                _in.write(';')
                _in.write(str(arr_y[j]))
                _in.write('\n')


if __name__ == '__main__':
    main()
