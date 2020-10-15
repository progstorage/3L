//#include "stdafx.h"
#include <stdlib.h>
#include <cmath> 
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;
using vector_pairs = vector<pair<double, double>>;
using Matrix = vector<vector<double>>;


vector_pairs sort(const vector<double>& x, const vector<double>& y) {
	//функция сортирует исходные координаты в смысле евклидовой метрики
	// и возвращает вектор пар

	int len = x.size();
	vector_pairs sorted_points(len);
	sorted_points[0].first = x[0];
	sorted_points[0].second = y[0];
	vector<int> flags(len);

	/*double **sorted_point = new double*[len];
	for (int i = 0; i < len; i++) {
		sorted_point[i] = new double[2];
	}*/
	//int* flag = new int[len];

	for (int i = 0; i < len; i++) { 
		flags[i] = 0;
	}
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
		if (flags[i] == 0) {
			for (int j = 0; j < 11; j++) {
				if ((j != i) && (flags[j] == 0)) {
					dist = sqrt((x[j] - tmp_x)*(x[j] - tmp_x) + (y[j] - tmp_y)*(y[j] - tmp_y));
					//cout << "j = " << j << ' ' << dist << endl;
					if (dist <= min_dist) {
						min_dist = dist;
						min_dist_ind = j;
					}
				}
				//else { continue; }
			}
			sorted_points[k + 1].first = x[min_dist_ind];
			sorted_points[k + 1].second = y[min_dist_ind];
			flags[i] = 1;
			i = min_dist_ind;
			min_dist = 100;
		}
		//else { continue; }
	}

	return sorted_points;
}


void levels(vector_pairs points) {
	// функция принимает на вход пары точек
	// и записывает полученные уровни в txt файл

	int len = points.size();
	vector_pairs l1(len);
	vector_pairs l3(len);

	//double **l1 = new double*[len];
	//double **l3 = new double*[len];
	/*for (int i = 0; i < len; i++) {
		l1[i] = new double[2];
		l3[i] = new double[2];
	}*/

	double tmpx, tmpy, n;
	double d = 0.08;
	for (int i = 1; i < len - 1; i++) {
		tmpx = points[i + 1].second - points[i - 1].second;
		tmpy = -(points[i + 1].first - points[i - 1].first);
		n = sqrt((tmpx*tmpx + tmpy * tmpy));
		tmpx = tmpx / n;
		tmpy = tmpy / n;
		l1[i].first = points[i].first + d * tmpx;
		l1[i].second = points[i].second + d * tmpy;
		l3[i].first = points[i].first - d * tmpx;
		l3[i].second = points[i].second - d * tmpy;
		if (i == 1) {
			l1[0].first = points[0].first + d * tmpx;
			l1[0].second = points[0].second + d * tmpy;
			l3[0].first = points[0].first - d * tmpx;
			l3[0].second = points[0].second - d * tmpy;
		}
		if (i == len - 2) {
			l1[len - 1].first = points[len - 1].first + d * tmpx;
			l1[len - 1].second = points[len - 1].second + d * tmpy;
			l3[len - 1].first = points[len - 1].first - d * tmpx;
			l3[len - 1].second = points[len - 1].second - d * tmpy;
		}
	}
	ofstream outfile ("out.txt");
	//fstream f;
	//f.open("out.txt", fstream::in | fstream::out);
	for (int i = 0; i < len; i++) {
		outfile << fixed << l1[i].first << "," << l1[i].second << "," << l3[i].first << "," << l3[i].second << endl;
	}
	outfile.close();
	/*for (int i = 0; i < 10; i++)
	{
		delete[] l1[i];
		delete[] l3[i];
	}
	delete[] l1;
	delete[] l3;*/
}


void print_Matrix(const Matrix& M) {
	// вывести матрицу
	for (const auto& v : M) {
		for (const auto& e : v)
			cout << e << " ";
		cout << endl;
	}
}


Matrix transpose_Matrix(const Matrix& M) {
	// транспонирование матрицы
	int n = M.size();
	int m = M[0].size();
	Matrix N(m, vector<double>(n));

	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
			N[j][i] = M[i][j];
	
	return N;
}


Matrix mult_Matrix(const Matrix& M, const Matrix& N) {
	// перемножение двух матриц
	int m1 = M.size();
	int m2 = M[0].size();
	int n1 = N.size();
	int n2 = N[0].size();
	if (m2 != n1) {
		cerr << "Wrong sizes!";
		exit(0);
	}
	Matrix K(m1, vector<double>(n2));

	double tmp;
	for (int i = 0; i < m1; i++) {
		for (int j = 0; j < n2; j++) {
			tmp = 0.0;
			//cout << i << " " << j << endl;
			for (int k = 0; k < m2; k++) {
				//cout << "k: " << k << endl;
				tmp += M[i][k] * N[k][j];
			}			
			K[i][j] = tmp;
		}
	}
	
	return K;
}


void fill_3L_Matrix_2nd_power(Matrix& M, const vector_pairs& points) {
	// заполняет матрицу (2я степень неявной функции)
	for (int i = 0; i < M.size(); i++) {
		M[i][0] = 1;
		M[i][1] = points[i].first;
		M[i][2] = points[i].second;
		M[i][3] = points[i].first * points[i].first;
		M[i][4] = points[i].first * points[i].second;
		M[i][5] = points[i].second * points[i].second;
	}
}


double det_Matrix(const Matrix& M) {
	// находит определитель квадратной матрицы
}


Matrix inverse_Matrix(const Matrix& M) {
	// возвращает обратную матрицу
}


int main(void) {
	//            [-1, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0] 
	//            [ 0,  0.6,  0.8,  0.9,  1.0, 1.0, 1.0, 0.9, 0.8, 0.6, 0.0]

	// test data
	vector<double> x = { -1, -0.8, -0.6, -0.4, 0.8, -0.2, 0.0, 0.2, 1.0, 0.4, 0.6 };
	vector<double> y = {  0,  0.6,  0.8,  0.9, 0.6,  1.0, 1.0, 1.0, 0.0, 0.9, 0.8 };

	int len = x.size();
	vector_pairs sorted_point = sort(x, y);

	levels(sorted_point);
	system("python plot_3L.py");		// отрисовка
	
	Matrix M(len, vector<double>(6));
	fill_3L_Matrix_2nd_power(M, sorted_point);

	sorted_point.clear();

	return(0);
}
