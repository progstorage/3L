#pragma once
//#include "stdafx.h"

#include <stdlib.h>
#include <cmath> 
#include <iostream>
#include <fstream>
#include <filesystem>
#include <string>
#include <vector>
#include "windows.h"
#include <ctime>
#include <omp.h>


// 1 - тестовый режим; 0 - обычный
#define TEST_MODE 1

using namespace std;
namespace fs = experimental::filesystem;

using vector_pairs =	vector<pair<double, double>>;
using Matrix =			vector<vector<double>>;

vector_pairs	sort(const vector<double>&, const vector<double>&, double*);
Matrix			transpose_Matrix(const Matrix&);
Matrix			mult_Matrix(const Matrix&, const Matrix&);
Matrix			mult_Matrix_multithread(const Matrix&, const Matrix&);
Matrix			inverse_Matrix(const Matrix&);
Matrix			generate_random_Matrix(const int, const int);
void			fill_3L_Matrix_2nd_power(Matrix&, const vector_pairs&, const vector_pairs&, const vector_pairs&);
void			levels(vector_pairs, double, vector_pairs&, vector_pairs&);
void			print_Matrix(const Matrix&);
void			solve_system(Matrix&, double);
double			det_Matrix(const Matrix&);
double			diff(pair<double, double>, pair<double, double>);


class Menu {
public:
	int execution_method; //1 - старый метод соединения точек, 2 - новый метод с разделением точек на группы
	string dir;
	vector<string> files;
	int n;
	vector<double> x;
	vector<double> y;
	string test_file;
	bool to_draw;

	Menu(string);
	void get_files(string);
	void print_menu();
	void select_file();
	void split(string, double*, double*);
	void read_file();
	void start();
	void select_execution_method();
};
