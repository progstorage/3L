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
vector<vector_pairs>	points_clusters_array;		//вектор групп точек
int						points_clusters_count = 0;	//количество групп точек


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
class cubic_spline
{
private:
	// Структура, описывающая сплайн на каждом сегменте сетки
	struct spline_tuple
	{
		double a, b, c, d, x;
	};
	spline_tuple *splines; // Сплайн
	std::size_t n; // Количество узлов сетки
	void free_mem(); // Освобождение памяти
public:
	cubic_spline(); //конструктор
	~cubic_spline(); //деструктор
	// Построение сплайна
	// x - узлы сетки, должны быть упорядочены по возрастанию, кратные узлы запрещены
	// y - значения функции в узлах сетки
	// n - количество узлов сетки
	void build_spline(vector<double>& x, vector<double>& y, std::size_t n);
	void write_poinst(cubic_spline Spline, double min, double max, int numpoints); // Записывает значения сплайна в точках промежутка [xmin-3, xmax+3]
	double f(double x) const;	// Вычисление значения интерполированной функции в произвольной точке
};