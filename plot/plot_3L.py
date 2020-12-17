# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
import os


def plot_g():
    x = []
    y = []
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


    l1 = [] #внутренний уровень
    l3 = [] #внешний уровень
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

    a = []
    for file in os.listdir("coefs/"):
        tmpa = []
        with open("coefs/" + file, 'r') as f:
            for line in f:
                tmpa.append(float(line))
        a.append(tmpa)


    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.grid(True, which='both')
    ax.set_ylim(min_y-1, max_y+1)
    ax.set_xlim(min_x-1, max_x+1)

    for i in range(len(x)):
        ax.scatter(x[i], y[i])
        ax.plot(x[i], y[i])
    
    for j in range(len(l1)):
        ax.scatter([i[0] for i in l1[j]], [i[1] for i in l1[j]])
        ax.plot([i[0] for i in l1[j]], [i[1] for i in l1[j]])

    for j in range(len(l3)):
        plt.scatter([i[0] for i in l3[j]], [i[1] for i in l3[j]])
        plt.plot([i[0] for i in l3[j]], [i[1] for i in l3[j]])

    for j in range(len(a)):
        f = lambda x, y: a[j][0]+a[j][1]*x+a[j][2]*y+a[j][3]*x*x+a[j][4]*x*y+a[j][5]*y*y
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