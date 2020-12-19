#include <vector>

class Matrix {
private:
	std::vector<std::vector<double> > elem;
public:
	Matrix(int, int);
	Matrix(const Matrix&);
	~Matrix();
	void resize(int, int);
	int rows_num() const;
	int cols_num() const;
	double det();
	void print();
	void fill_3L_Matrix_2nd_power(const vector_pairs&, const vector_pairs&, const vector_pairs&);
	Matrix transpose();
	Matrix inverse();
	Matrix& operator*(const Matrix&) const; // многопоточное умножение
	double* operator[](int) const; // многопоточное умножение
};