// ConsoleApplication1.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"
#include <cmath> 
#include <iostream>
#include <fstream>

using namespace std;

double** sort(double* x, double* y, int len) {
	double **sorted_point = new double*[len];
	for (int i = 0; i < len; i++) {
		sorted_point[i] = new double[2];
	}
	sorted_point[0][0] = x[0];
	sorted_point[0][1] = y[0];
	int* flag = new int[len];
	for (int i = 0; i < len; i++) { flag[i] = 0; }
	double dist = 0;
	double min_dist = sqrt((x[1] - x[0])*(x[1] - x[0]) + (y[1] - y[0])*(y[1] - y[0]));
	int min_dist_ind = 0;
	int i = 0;
	double tmp_x = 0;
	double tmp_y = 0;
	for (int k = 0; k < len - 1; k++) {
		//cout << "i = " << i << endl;
		tmp_x = x[i];
		tmp_y = y[i];
		if (flag[i] == 0) {
			for (int j = 0; j < 11; j++) {
				if ((j != i) && (flag[j] == 0)) {
					dist = sqrt((x[j] - tmp_x)*(x[j] - tmp_x) + (y[j] - tmp_y)*(y[j] - tmp_y));
					//cout << "j = " << j << ' ' << dist << endl;
					if (dist <= min_dist) {
						min_dist = dist;
						min_dist_ind = j;
					}
				}
				else { continue; }
			}
			sorted_point[k + 1][0] = x[min_dist_ind];
			sorted_point[k + 1][1] = y[min_dist_ind];
			flag[i] = 1;
			i = min_dist_ind;
			min_dist = 100;
		}
		else { continue; }
	}

	return sorted_point;
}

void levels(double** points, int len) {
	double **l1 = new double*[len];
	double **l3 = new double*[len];
	for (int i = 0; i < len; i++) {
		l1[i] = new double[2];
		l3[i] = new double[2];
	}

	double tmpx, tmpy, n;
	double d = 0.08;
	for (int i = 1; i < len - 1; i++) {
		tmpx = points[i + 1][1] - points[i - 1][1];
		tmpy = -(points[i + 1][0] - points[i - 1][0]);
		n = sqrt((tmpx*tmpx + tmpy * tmpy));
		tmpx = tmpx / n;
		tmpy = tmpy / n;
		l1[i][0] = points[i][0] + d * tmpx;
		l1[i][1] = points[i][1] + d * tmpy;
		l3[i][0] = points[i][0] - d * tmpx;
		l3[i][1] = points[i][1] - d * tmpy;
		if (i == 1) {
			l1[0][0] = points[0][0] + d * tmpx;
			l1[0][1] = points[0][1] + d * tmpy;
			l3[0][0] = points[0][0] - d * tmpx;
			l3[0][1] = points[0][1] - d * tmpy;
		}
		if (i == len - 2) {
			l1[len - 1][0] = points[len - 1][0] + d * tmpx;
			l1[len - 1][1] = points[len - 1][1] + d * tmpy;
			l3[len - 1][0] = points[len - 1][0] - d * tmpx;
			l3[len - 1][1] = points[len - 1][1] - d * tmpy;
		}
	}
	fstream f;
	f.open("out.txt", fstream::in | fstream::out);
	for (int i = 0; i < len; i++) {
		f << fixed << l1[i][0] << "," << l1[i][1] << "," << l3[i][0] << "," << l3[i][1] << endl;
	}
	f.close();
	for (int i = 0; i < 10; i++)
	{
		delete[] l1[i];
		delete[] l3[i];
	}
	delete[] l1;
	delete[] l3;
}

int main() {
	//            [-1, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0] 
	//            [ 0,  0.6,  0.8,  0.9,  1.0, 1.0, 1.0, 0.9, 0.8, 0.6, 0.0]

	double x[] = { -1, -0.8, -0.6, -0.4, 0.8, -0.2, 0.0, 0.2, 1.0, 0.4, 0.6 };
	double y[] = { 0,  0.6,  0.8,  0.9, 0.6,  1.0, 1.0, 1.0, 0.0, 0.9, 0.8 };

	int len = sizeof(x) / sizeof(x[0]);
	double** sorted_point = sort(x, y, len);

	levels(sorted_point, len);
	system("python f1.py");

	for (int i = 0; i < 10; i++)
	{
		delete[] sorted_point[i];
	}
	delete[] sorted_point;
	return(0);
}
