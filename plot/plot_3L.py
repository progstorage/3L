# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
import os


def plot_g():
    x = []
    y = []
    try:
        for file in os.listdir("sorted_points/"):
            tmpx = []
            tmpy = []
            with open("sorted_points/" + file, "r") as f:
                for line in f:
                    l = line.strip('\n').split(',')
                    tmpx.append(float(l[0]))
                    tmpy.append(float(l[1]))
            x.append(tmpx)
            y.append(tmpy)
        minx = []
        maxx = []
        for i in x:
            minx.append(min(i))
            maxx.append(max(i))
        min_x = min(minx)
        max_x = max(maxx)
        miny = []
        maxy = []
        for i in y:
            miny.append(min(i))
            maxy.append(max(i))
        min_y = min(miny)
        max_y = max(maxy)
    except Exception as exc:
        print(f"Can't plot, error: {exc}")
        return

    l1 = [] #внутренний уровень
    l3 = [] #внешний уровень
    try:
        for file in os.listdir("levels/"):
            tmpl1 = []
            tmpl3 = []
            with open("levels/" + file, 'r') as f:
                for line in f:
                    l = line.strip('\n').split(',')
                    tmpl1.append([float(l[0]), float(l[1])])
                    tmpl3.append([float(l[2]), float(l[3])])
            l1.append(tmpl1)
            l3.append(tmpl3)
    except Exception as exc:
        print(f"Can't plot levels, error: {exc}")

    a = []
    try:
        for file in os.listdir("coefs/"):
            with open("coefs/" + file, 'r') as f:
                for line in f:
                    a.append(float(line))
    except Exception as exc:
        print(f"Can't plot polinom, error: {exc}")

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.grid(True, which='both')
    ax.set_ylim(min_y-1, max_y+1)
    ax.set_xlim(min_x-1, max_x+1)

    if x:
        ax.scatter(x, y)
    if l1:
        for j in range(len(l1)):
            ax.scatter([i[0] for i in l1[j]], [i[1] for i in l1[j]], color="green")
            plt.scatter([i[0] for i in l3[j]], [i[1] for i in l3[j]], color="orange")

    if a:
        len_a = len(a)
        if len_a == 6:
            f = lambda x, y: a[0]+a[1]*x+a[2]*y+a[3]*x*x+a[4]*x*y+a[5]*y*y
        elif len_a == 10:
            f = lambda x, y: a[0]+a[1]*x+a[2]*y+a[3]*x*x+a[4]*x*y+a[5]*y*y +a[6]*x**3+a[7]*x**2*y+a[8]*x*y**2+a[9]*y**3
        elif len_a == 15:
            f = lambda x, y: a[0]+a[1]*x+a[2]*y+a[3]*x*x+a[4]*x*y+a[5]*y*y +a[6]*x**3+a[7]*x**2*y+a[8]*x*y**2+a[9]*y**3+a[10]*x**4+a[11]*x**3*y+a[12]*x**2*y**2+a[13]*x*y**3+a[14]*y**4
        elif len_a == 21:
            f = lambda x, y: a[0]+a[1]*x+a[2]*y+a[3]*x*x+a[4]*x*y+a[5]*y*y +a[6]*x**3+a[7]*x**2*y+a[8]*x*y**2+a[9]*y**3+a[10]*x**4+a[11]*x**3*y+a[12]*x**2*y**2+a[13]*x*y**3+a[14]*y**4+a[15]*x**5+a[16]*x**4*y+a[17]*x**3*y**2+a[18]*x**2*y**3+a[19]*x*y**4+a[20]*y**5
        elif len_a == 28:
            f = lambda x, y: a[0]+a[1]*x+a[2]*y+a[3]*x*x+a[4]*x*y+a[5]*y*y +a[6]*x**3+a[7]*x**2*y+a[8]*x*y**2+a[9]*y**3+a[10]*x**4+a[11]*x**3*y+a[12]*x**2*y**2+a[13]*x*y**3+a[14]*y**4+a[15]*x**5+a[16]*x**4*y+a[17]*x**3*y**2+a[18]*x**2*y**3+a[19]*x*y**4+a[20]*y**5+a[21]*x**6+a[22]*x**5*y+a[23]*x**4*y**2+a[24]*x**3*y**3+a[25]*x**2*y**4+a[26]*x*y**5+a[27]*y**6
        elif len_a == 36:
            f = lambda x, y: a[0]+a[1]*x+a[2]*y+a[3]*x*x+a[4]*x*y+a[5]*y*y +a[6]*x**3+a[7]*x**2*y+a[8]*x*y**2+a[9]*y**3+a[10]*x**4+a[11]*x**3*y+a[12]*x**2*y**2+a[13]*x*y**3+a[14]*y**4+a[15]*x**5+a[16]*x**4*y+a[17]*x**3*y**2+a[18]*x**2*y**3+a[19]*x*y**4+a[20]*y**5+a[21]*x**6+a[22]*x**5*y+a[23]*x**4*y**2+a[24]*x**3*y**3+a[25]*x**2*y**4+a[26]*x*y**5+a[27]*y**6+a[28]*x**7+a[29]*x**6*y+a[30]*x**5*y**2+a[31]*x**4*y**3+a[32]*x**3*y**4+a[33]*x**2*y**5+a[34]*x*y**6+a[35]*y**7
        delta = 0.025
        x_range = np.arange(min_y-1, max_y+1, delta)
        y_range = np.arange(min_x-1, max_x+1, delta)
        X, Y = np.meshgrid(x_range, y_range)
        F = f(X, Y)
        ax.contour(X, Y, F, [0], colors = ['red'])
        #plt.gca().set_aspect('equal', adjustable='box')
        #plt.draw()
    plt.show()
    

if __name__ == '__main__':
    plot_g();