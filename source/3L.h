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

using namespace std;
namespace fs = experimental::filesystem;
using vector_pairs = vector<pair<double, double>>;
using Matrix = vector<vector<double>>;


double diff(pair<double, double>, pair<double, double>);
vector_pairs sort(const vector<double>&, const vector<double>&, double*);
void levels(vector_pairs, double, vector_pairs&, vector_pairs&);
void print_Matrix(const Matrix&);
Matrix transpose_Matrix(const Matrix&);
Matrix transpose_Matrix(const Matrix&);
Matrix mult_Matrix(const Matrix&, const Matrix&);
void fill_3L_Matrix_2nd_power(Matrix&, const vector_pairs&, const vector_pairs&, const vector_pairs&);
double det_Matrix(const Matrix&);
Matrix inverse_Matrix(const Matrix&);
void solve_system(Matrix&, double);


class Menu {
public:
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
};
