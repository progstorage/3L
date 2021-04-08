#ifndef MATRIXX
#define MATRIXX

#include <vector>
#include <iostream>
#include <cmath>
using namespace std;
using vector_pairs = std::vector<std::pair<double, double> >;


class Matrix {
private:
	std::vector<vector<double>> elem;
public:
	Matrix() = default;
	Matrix(int, int);
	Matrix(const Matrix&);
	void resize(int, int);
	int rows_num() const;
	int cols_num() const;
	double det();
	void print();
	void fill_3L_Matrix_2nd_power(const vector_pairs&, const vector_pairs&, const vector_pairs&);
	void fill_3L_Matrix_3d_power(const vector_pairs&, const vector_pairs&, const vector_pairs&);
	void fill_3L_Matrix_4th_power(const vector_pairs&, const vector_pairs&, const vector_pairs&);
	void fill_3L_Matrix_5th_power(const vector_pairs&, const vector_pairs&, const vector_pairs&);
	void fill_3L_Matrix_6th_power(const vector_pairs&, const vector_pairs&, const vector_pairs&);
	void fill_3L_Matrix_7th_power(const vector_pairs&, const vector_pairs&, const vector_pairs&);
	Matrix transpose() const;	// returns transposed copy of matrix
	Matrix inverse();
	Matrix element_power(double);
	Matrix element_div(Matrix&);
	Matrix element_div(double);
	
	
	vector<double>& operator[](int);
	const vector<double>& operator[](int) const;
	Matrix& operator=(const Matrix&);

	//double& operator()(int, int);
	friend Matrix mult_Matrix(const Matrix&, const Matrix&);	// однопоточное умножение
	friend Matrix operator*(const Matrix&, const Matrix&);		// многопоточное умножение
	friend Matrix operator*(const double, const Matrix&);
	friend bool operator==(const Matrix&, const Matrix&);
	friend Matrix operator-(const Matrix&);
	friend Matrix operator-(const Matrix&, const Matrix&);
	friend Matrix operator+(const Matrix&, const Matrix&);
};

#endif