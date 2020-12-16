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
#include "config.h"

using namespace std;

#ifdef __APPLE__
	namespace fs = __fs::filesystem;
#else
	namespace fs = experimental::filesystem;
#endif

using vector_pairs =	vector<pair<double, double>>;
using Matrix =			vector<vector<double>>;

Matrix		transpose_Matrix(const Matrix&);
Matrix		mult_Matrix(const Matrix&, const Matrix&);
Matrix		mult_Matrix_multithread(const Matrix&, const Matrix&);
Matrix		inverse_Matrix(const Matrix&);
Matrix		generate_random_Matrix(const int, const int);
void		fill_3L_Matrix_2nd_power(Matrix&, const vector_pairs&, const vector_pairs&, const vector_pairs&);
void		levels(vector_pairs, double, vector_pairs&, vector_pairs&, int);
void		print_Matrix(const Matrix&);
void		solve_system(Matrix&, double, int);
void		fill_levels_multithread(double);
double		det_Matrix(const Matrix&);
double		diff(pair<double, double>, pair<double, double>);
int			sort_method_2(const vector<double>&, const vector<double>&, double*);
int			old_sort(const vector<double>&, const vector<double>&, double*);

// global vars
vector<vector_pairs>	points_clusters_array;		//������ ����� �����
int						points_clusters_count = 0;	//���������� ����� �����


class Menu {
public:
	int execution_method; //1 - ������ ����� ���������� �����, 2 - ����� ����� � ����������� ����� �� ������
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
class cubic_spline
{
private:
	// ���������, ����������� ������ �� ������ �������� �����
	struct spline_tuple
	{
		double a, b, c, d, x;
	};
	spline_tuple *splines; // ������
	std::size_t n; // ���������� ����� �����
	void free_mem(); // ������������ ������
public:
	cubic_spline(); //�����������
	~cubic_spline(); //����������
	// ���������� �������
	// x - ���� �����, ������ ���� ����������� �� �����������, ������� ���� ���������
	// y - �������� ������� � ����� �����
	// n - ���������� ����� �����
	void build_spline(vector<double>& x, vector<double>& y, std::size_t n);
	void write_poinst(cubic_spline Spline, double min, double max, int numpoints); // ���������� �������� ������� � ������ ���������� [xmin-3, xmax+3]
	double f(double x) const;	// ���������� �������� ����������������� ������� � ������������ �����
};