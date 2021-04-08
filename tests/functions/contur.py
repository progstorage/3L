import numpy as np
import matplotlib.pyplot as plt
import os
import cv2

num = 20

img = cv2.imread('contur_pictures/' + str(num) + '.TIFF')
edge = cv2.Canny(img, 100, 200)
arr_x = []
arr_y = []
for y in range(0, edge.shape[0], 7):
    for x in range(0, edge.shape[1], 7):
        if edge[y, x] != 0:
            arr_x.append(x)
            arr_y.append(y)

file_name = '../in_2/in_070.txt'
with open(file_name, 'w+') as _in:
    for j in range(len(arr_x)):
        _in.write(str(arr_x[j]))
        _in.write(';')
        _in.write(str(arr_y[j]))
        _in.write('\n')

