#pragma once
//#include "stdafx.h"

#ifdef _WIN32
	#include "windows.h"
	#include <process.h>
#endif

#include <stdlib.h>
#include <cmath> 
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <limits>
#include <filesystem>
// #include "source/config.h"

using namespace std;
using vector_pairs = std::vector<std::pair<double, double> >;

#include "Matrix.cpp"
#include "cubic_spline.cpp"

#ifdef __APPLE__
	namespace fs = __fs::filesystem;
#else
	namespace fs = experimental::filesystem;
#endif

Matrix		generate_random_Matrix(const int, const int);
void		levels(vector_pairs, double, vector_pairs&, vector_pairs&, int);
void		solve_system(Matrix&, double, int);
void		fill_levels_multithread(double);
double		diff(pair<double, double>, pair<double, double>);
int			sort_method_2(const vector<double>&, const vector<double>&, double*);
int			old_sort(const vector<double>&, const vector<double>&, double*);

// global vars
vector<vector_pairs>	points_clusters_array;		//вектор групп точек
int						points_clusters_count = 0;	//количество групп точек