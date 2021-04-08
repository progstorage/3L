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
//#include <ctime>
#include <cstdlib>
#include <limits>
#include <chrono>

using namespace std;
using vector_pairs = vector<pair<double, double> >;
#define _SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING
#ifdef __APPLE__
#include <filesystem>
namespace fs = __fs::filesystem;
#else
#include <experimental/filesystem>
namespace fs = experimental::filesystem;
#endif

#include "Matrix.cpp"
#include "cubic_spline.cpp"

Matrix		generate_random_Matrix(const int, const int);
void		levels(vector_pairs, double, vector_pairs&, vector_pairs&, int);
void		test_spline_norm_levels(vector_pairs, double, vector_pairs&, vector_pairs&, int);
void		solve_system(Matrix&, double, int);
void		fill_levels_multithread(double);
double		diff(pair<double, double>, pair<double, double>);
int			sort_method_2(const vector<double>&, const vector<double>&, double*);
int			old_sort(const vector<double>&, const vector<double>&, double*);

// new alg
vector<double>	matr_max(const Matrix&);
int				vec_min_ind(vector<double>);
double			vec_max(vector<double>);
double			vec_min(vector<double>);
double			vec_sum(vector<double>, int);
vector<double>	vec_mult(vector<double>, vector<double>);
double			mean(vector<double>);
void			quickSort(vector<double>&, int, int, vector<int>&);
void			mls(double*, double*, double*, vector<double>, vector<double>);
double			pre3L(vector<double>&, vector<double>&);


// global vars
vector<vector_pairs>	points_clusters_array;		//âåêòîð ãðóïï òî÷åê
int						points_clusters_count = 0;	//êîëè÷åñòâî ãðóïï òî÷åê