#ifndef MATRIXX
#define MATRIXX

#include <vector>
#include <iostream>

using namespace std;
using vector_pairs = std::vector<std::pair<double, double> >;


class Matrix {
private:
	std::vector<vector<double>> elem;
public:
	Matrix(int, int);
	Matrix(Matrix&);
	~Matrix();
	void resize(int, int);
	int rows_num() const;
	int cols_num() const;
	double det();
	void print();
	void fill_3L_Matrix_2nd_power(const vector_pairs&, const vector_pairs&, const vector_pairs&);
	Matrix transpose();
	Matrix inverse();
	Matrix& operator*(Matrix&) const; // многопоточное умножение
	
	vector<double>& operator[](int);

	//double& operator()(int, int);
	friend Matrix mult_Matrix(const Matrix&, const Matrix&);
};

#endif