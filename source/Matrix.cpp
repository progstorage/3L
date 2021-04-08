#include "Matrix.h"
// конструктор
Matrix::Matrix(int rows, int cols) {
	this->elem.resize(rows, vector<double>(cols));
}

Matrix::Matrix(const Matrix& M) {
	int rows = M.rows_num();
	int cols = M.cols_num();
	this->elem.resize(rows, vector<double>(cols));
	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
			elem[i][j] = M[i][j];
}

// возврат размерности
int Matrix::rows_num() const { return elem.size(); }
int Matrix::cols_num() const { return elem[0].size(); }

// поиск определител€
double Matrix::det() {
	// определитель ћетод √аусса
	const double EPS = 1E-9;
	int n = elem.size();
	Matrix N(*this);	//  опирование матрицы
	double det = 1;
	for (int i = 0; i < n; ++i) {
		int k = i;
		for (int j = i + 1; j < n; ++j) {
			if (abs(N[j][i]) > abs(N[k][i]))
				k = j;
		}
		if (abs(N[k][i]) < EPS) {
			det = 0;
			return det;
		}
		swap(N[i], N[k]);
		for (int j = 0; j < n; ++j)
			if (j != i && abs(N[j][i]) > EPS)
				for (int k = i + 1; k < n; ++k)
					N[j][k] -= N[i][k] * N[j][i];
	}
	return det;
}

// печать матрицы в консоль
void Matrix::print() {
	for (const auto& v : elem) {
		for (const auto& e : v)
			cout << e << "\t";
		cout << endl;
	}
}

bool operator==(const Matrix& M, const Matrix& N) {
	int m1 = M.rows_num();
	int m2 = M.cols_num();
	int n1 = N.rows_num();
	int n2 = N.cols_num();
	if (m1 != n1 or m2 != n2) {
		cerr << "Wrong sizes for comparison!";
		exit(0);
	}

	for (int i = 0; i < m1; i++)
		for (int j = 0; j < m2; j++)
			if (M[i][j] != N[i][j])
				return false;
	return true;
}

void Matrix::fill_3L_Matrix_power(const vector_pairs& points, const vector_pairs& l1, const vector_pairs& l3, int power) {
	int max_col = 0;
	for (int i = 1; i <= power+1; i++) {
		max_col += i;
	}
	int cur_row; int max_pow;
	for (int i = 0; i < this->elem.size() / 3; i++) {
		cur_row = 0; max_pow = 1;
		while (cur_row < max_col) {
			for (int b = 0; b < max_pow; b++) {
				this->elem[i][cur_row] = pow(l1[i].first, max_pow - 1 - b) * pow(l1[i].second, 0 + b);
				cur_row++;
			}
			max_pow++;
		}
	}
	int j = 0;
	for (int i = this->elem.size() / 3; i < 2 * this->elem.size() / 3; i++) {
		cur_row = 0; max_pow = 1;
		while (cur_row < max_col) {
			for (int b = 0; b < max_pow; b++) {
				this->elem[i][cur_row] = pow(points[j].first, max_pow - 1 - b) * pow(points[j].second, 0 + b);
				cur_row++;
			}
			max_pow++;
		}
		j++;
	}
	j = 0;
	for (int i = 2 * this->elem.size() / 3; i < this->elem.size(); i++) {
		cur_row = 0; max_pow = 1;
		while (cur_row < max_col) {
			for (int b = 0; b < max_pow; b++) {
				this->elem[i][cur_row] = pow(l3[j].first, max_pow - 1 - b) * pow(l3[j].second, 0 + b);
				cur_row++;
			}
			max_pow++;
		}
		j++;
	}
}

// заполнение матрицы
void Matrix::fill_3L_Matrix_2nd_power(const vector_pairs& points, const vector_pairs& l1, const vector_pairs& l3) {
	// заполн€ет матрицу (2€ степень не€вной функции)
	// ћатрица вида 
	// [M-]
	// [M0]
	// [M+]
	for (int i = 0; i < this->elem.size() / 3; i++) {
		this->elem[i][0] = 1;
		this->elem[i][1] = l1[i].first;
		this->elem[i][2] = l1[i].second;
		this->elem[i][3] = l1[i].first * l1[i].first;
		this->elem[i][4] = l1[i].first * l1[i].second;
		this->elem[i][5] = l1[i].second * l1[i].second;
	}
	int j = 0;
	for (int i = this->elem.size() / 3; i < 2 * this->elem.size() / 3; i++) {
		this->elem[i][0] = 1;
		this->elem[i][1] = points[j].first;
		this->elem[i][2] = points[j].second;
		this->elem[i][3] = points[j].first * points[j].first;
		this->elem[i][4] = points[j].first * points[j].second;
		this->elem[i][5] = points[j].second * points[j].second;
		j++;
	}
	j = 0;
	for (int i = 2 * this->elem.size() / 3; i < this->elem.size(); i++) {
		this->elem[i][0] = 1;
		this->elem[i][1] = l3[j].first;
		this->elem[i][2] = l3[j].second;
		this->elem[i][3] = l3[j].first * l3[j].first;
		this->elem[i][4] = l3[j].first * l3[j].second;
		this->elem[i][5] = l3[j].second * l3[j].second;
		j++;
	}
}


void Matrix::fill_3L_Matrix_3d_power(const vector_pairs& points, const vector_pairs& l1, const vector_pairs& l3) {
	// заполн€ет матрицу (3€ степень не€вной функции)
	// ћатрица вида 
	// [M-]
	// [M0]
	// [M+]
	for (int i = 0; i < this->elem.size() / 3; i++) {
		this->elem[i][0] = 1;
		this->elem[i][1] = l1[i].first;
		this->elem[i][2] = l1[i].second;
		this->elem[i][3] = l1[i].first * l1[i].first;
		this->elem[i][4] = l1[i].first * l1[i].second;
		this->elem[i][5] = l1[i].second * l1[i].second;

		this->elem[i][6] = l1[i].first * l1[i].first * l1[i].first;
		this->elem[i][7] = l1[i].first * l1[i].first * l1[i].second;
		this->elem[i][8] = l1[i].second * l1[i].second * l1[i].first;
		this->elem[i][9] = l1[i].second * l1[i].second * l1[i].second;
	}
	int j = 0;
	for (int i = this->elem.size() / 3; i < 2 * this->elem.size() / 3; i++) {
		this->elem[i][0] = 1;
		this->elem[i][1] = points[j].first;
		this->elem[i][2] = points[j].second;
		this->elem[i][3] = points[j].first * points[j].first;
		this->elem[i][4] = points[j].first * points[j].second;
		this->elem[i][5] = points[j].second * points[j].second;

		this->elem[j][6] = points[j].first * points[j].first * points[j].first;
		this->elem[j][7] = points[j].first * points[j].first * points[j].second;
		this->elem[j][8] = points[j].second * points[j].second * points[j].first;
		this->elem[j][9] = points[j].second * points[j].second * points[j].second;
		j++;
	}
	j = 0;
	for (int i = 2 * this->elem.size() / 3; i < this->elem.size(); i++) {
		this->elem[i][0] = 1;
		this->elem[i][1] = l3[j].first;
		this->elem[i][2] = l3[j].second;
		this->elem[i][3] = l3[j].first * l3[j].first;
		this->elem[i][4] = l3[j].first * l3[j].second;
		this->elem[i][5] = l3[j].second * l3[j].second;

		this->elem[j][6] = l3[j].first * l3[j].first * l3[j].first;
		this->elem[j][7] = l3[j].first * l3[j].first * l3[j].second;
		this->elem[j][8] = l3[j].second * l3[j].second * l3[j].first;
		this->elem[j][9] = l3[j].second * l3[j].second * l3[j].second;
		j++;
	}
}

void Matrix::fill_3L_Matrix_4th_power(const vector_pairs& points, const vector_pairs& l1, const vector_pairs& l3) {
	// заполн€ет матрицу (4 степень не€вной функции)
	// ћатрица вида 
	// [M-]
	// [M0]
	// [M+]
	for (int i = 0; i < this->elem.size() / 3; i++) {
		this->elem[i][0] = 1;
		this->elem[i][1] = l1[i].first;
		this->elem[i][2] = l1[i].second;
		
		this->elem[i][3] = l1[i].first * l1[i].first;
		this->elem[i][4] = l1[i].first * l1[i].second;
		this->elem[i][5] = l1[i].second * l1[i].second;

		this->elem[i][6] = l1[i].first * l1[i].first * l1[i].first;
		this->elem[i][7] = l1[i].first * l1[i].first * l1[i].second;
		this->elem[i][8] = l1[i].second * l1[i].second * l1[i].first;
		this->elem[i][9] = l1[i].second * l1[i].second * l1[i].second;
		
		this->elem[i][10] = l1[i].first * l1[i].first * l1[i].first * l1[i].first;
		this->elem[i][11] = l1[i].first * l1[i].first * l1[i].first * l1[i].second;
		this->elem[i][12] = l1[i].first * l1[i].second * l1[i].second * l1[i].first;
		this->elem[i][13] = l1[i].first * l1[i].second * l1[i].second * l1[i].second;
		this->elem[i][14] = l1[i].second * l1[i].second * l1[i].second * l1[i].second;
	}
	int j = 0;
	for (int i = this->elem.size() / 3; i < 2 * this->elem.size() / 3; i++) {
		this->elem[i][0] = 1;
		this->elem[i][1] = points[j].first;
		this->elem[i][2] = points[j].second;
		this->elem[i][3] = points[j].first * points[j].first;
		this->elem[i][4] = points[j].first * points[j].second;
		this->elem[i][5] = points[j].second * points[j].second;

		this->elem[j][6] = points[j].first * points[j].first * points[j].first;
		this->elem[j][7] = points[j].first * points[j].first * points[j].second;
		this->elem[j][8] = points[j].second * points[j].second * points[j].first;
		this->elem[j][9] = points[j].second * points[j].second * points[j].second;
		
		this->elem[i][10] = points[j].first * points[j].first * points[j].first * points[j].first;
		this->elem[i][11] = points[j].first * points[j].first * points[j].first * points[j].second;
		this->elem[i][12] = points[j].first * points[j].second * points[j].second * points[j].first;
		this->elem[i][13] = points[j].first * points[j].second * points[j].second * points[j].second;
		this->elem[i][14] = points[j].second * points[j].second * points[j].second * points[j].second;
		j++;
	}
	j = 0;
	for (int i = 2 * this->elem.size() / 3; i < this->elem.size(); i++) {
		this->elem[i][0] = 1;
		this->elem[i][1] = l3[j].first;
		this->elem[i][2] = l3[j].second;
		this->elem[i][3] = l3[j].first * l3[j].first;
		this->elem[i][4] = l3[j].first * l3[j].second;
		this->elem[i][5] = l3[j].second * l3[j].second;

		this->elem[j][6] = l3[j].first * l3[j].first * l3[j].first;
		this->elem[j][7] = l3[j].first * l3[j].first * l3[j].second;
		this->elem[j][8] = l3[j].second * l3[j].second * l3[j].first;
		this->elem[j][9] = l3[j].second * l3[j].second * l3[j].second;
		
		this->elem[i][10] = l3[j].first * l3[j].first * l3[j].first * l3[j].first;
		this->elem[i][11] = l3[j].first * l3[j].first * l3[j].first * l3[j].second;
		this->elem[i][12] = l3[j].first * l3[j].second * l3[j].second * l3[j].first;
		this->elem[i][13] = l3[j].first * l3[j].second * l3[j].second * l3[j].second;
		this->elem[i][14] = l3[j].second * l3[j].second * l3[j].second * l3[j].second;
		j++;
	}
}

void Matrix::fill_3L_Matrix_5th_power(const vector_pairs& points, const vector_pairs& l1, const vector_pairs& l3) {
	// заполн€ет матрицу (5 степень не€вной функции)
	// ћатрица вида 
	// [M-]
	// [M0]
	// [M+]
	for (int i = 0; i < this->elem.size() / 3; i++) {
		this->elem[i][0] = 1;
		this->elem[i][1] = l1[i].first;
		this->elem[i][2] = l1[i].second;
		
		this->elem[i][3] = l1[i].first * l1[i].first;
		this->elem[i][4] = l1[i].first * l1[i].second;
		this->elem[i][5] = l1[i].second * l1[i].second;

		this->elem[i][6] = l1[i].first * l1[i].first * l1[i].first;
		this->elem[i][7] = l1[i].first * l1[i].first * l1[i].second;
		this->elem[i][8] = l1[i].second * l1[i].second * l1[i].first;
		this->elem[i][9] = l1[i].second * l1[i].second * l1[i].second;
		
		this->elem[i][10] = l1[i].first * l1[i].first * l1[i].first * l1[i].first;
		this->elem[i][11] = l1[i].first * l1[i].first * l1[i].first * l1[i].second;
		this->elem[i][12] = l1[i].first * l1[i].second * l1[i].second * l1[i].first;
		this->elem[i][13] = l1[i].first * l1[i].second * l1[i].second * l1[i].second;
		this->elem[i][14] = l1[i].second * l1[i].second * l1[i].second * l1[i].second;
		
		this->elem[i][15] = l1[i].first * l1[i].first * l1[i].first * l1[i].first * l1[i].first;
		this->elem[i][16] = l1[i].first * l1[i].first * l1[i].first * l1[i].first * l1[i].second;
		this->elem[i][17] = l1[i].first * l1[i].first * l1[i].second * l1[i].second * l1[i].first;
		this->elem[i][18] = l1[i].first * l1[i].first * l1[i].second * l1[i].second * l1[i].second;
		this->elem[i][19] = l1[i].first * l1[i].second * l1[i].second * l1[i].second * l1[i].second;
		this->elem[i][20] = l1[i].second * l1[i].second * l1[i].second * l1[i].second * l1[i].second;
	}
	int j = 0;
	for (int i = this->elem.size() / 3; i < 2 * this->elem.size() / 3; i++) {
		this->elem[i][0] = 1;
		this->elem[i][1] = points[j].first;
		this->elem[i][2] = points[j].second;
		this->elem[i][3] = points[j].first * points[j].first;
		this->elem[i][4] = points[j].first * points[j].second;
		this->elem[i][5] = points[j].second * points[j].second;

		this->elem[j][6] = points[j].first * points[j].first * points[j].first;
		this->elem[j][7] = points[j].first * points[j].first * points[j].second;
		this->elem[j][8] = points[j].second * points[j].second * points[j].first;
		this->elem[j][9] = points[j].second * points[j].second * points[j].second;
		
		this->elem[i][10] = points[j].first * points[j].first * points[j].first * points[j].first;
		this->elem[i][11] = points[j].first * points[j].first * points[j].first * points[j].second;
		this->elem[i][12] = points[j].first * points[j].second * points[j].second * points[j].first;
		this->elem[i][13] = points[j].first * points[j].second * points[j].second * points[j].second;
		this->elem[i][14] = points[j].second * points[j].second * points[j].second * points[j].second;
		
		this->elem[i][15] = points[j].first * points[j].first * points[j].first * points[j].first * points[j].first;
		this->elem[i][16] = points[j].first * points[j].first * points[j].first * points[j].first * points[j].second;
		this->elem[i][17] = points[j].first * points[j].first * points[j].second * points[j].second * points[j].first;
		this->elem[i][18] = points[j].first * points[j].first * points[j].second * points[j].second * points[j].second;
		this->elem[i][19] = points[j].first * points[j].second * points[j].second * points[j].second * points[j].second;
		this->elem[i][20] = points[j].second * points[j].second * points[j].second * points[j].second * points[j].second;
		j++;
	}
	j = 0;
	for (int i = 2 * this->elem.size() / 3; i < this->elem.size(); i++) {
		this->elem[i][0] = 1;
		this->elem[i][1] = l3[j].first;
		this->elem[i][2] = l3[j].second;
		this->elem[i][3] = l3[j].first * l3[j].first;
		this->elem[i][4] = l3[j].first * l3[j].second;
		this->elem[i][5] = l3[j].second * l3[j].second;

		this->elem[j][6] = l3[j].first * l3[j].first * l3[j].first;
		this->elem[j][7] = l3[j].first * l3[j].first * l3[j].second;
		this->elem[j][8] = l3[j].second * l3[j].second * l3[j].first;
		this->elem[j][9] = l3[j].second * l3[j].second * l3[j].second;
		
		this->elem[i][10] = l3[j].first * l3[j].first * l3[j].first * l3[j].first;
		this->elem[i][11] = l3[j].first * l3[j].first * l3[j].first * l3[j].second;
		this->elem[i][12] = l3[j].first * l3[j].second * l3[j].second * l3[j].first;
		this->elem[i][13] = l3[j].first * l3[j].second * l3[j].second * l3[j].second;
		this->elem[i][14] = l3[j].second * l3[j].second * l3[j].second * l3[j].second;
		
		this->elem[i][15] = l3[j].first * l3[j].first * l3[j].first * l3[j].first * l3[j].first;
		this->elem[i][16] = l3[j].first * l3[j].first * l3[j].first * l3[j].first * l3[j].second;
		this->elem[i][17] = l3[j].first * l3[j].first * l3[j].second * l3[j].second * l3[j].first;
		this->elem[i][18] = l3[j].first * l3[j].first * l3[j].second * l3[j].second * l3[j].second;
		this->elem[i][19] = l3[j].first * l3[j].second * l3[j].second * l3[j].second * l3[j].second;
		this->elem[i][20] = l3[j].second * l3[j].second * l3[j].second * l3[j].second * l3[j].second;
		j++;
	}
}

void Matrix::fill_3L_Matrix_6th_power(const vector_pairs& points, const vector_pairs& l1, const vector_pairs& l3) {
	// заполн€ет матрицу (6 степень не€вной функции)
	// ћатрица вида 
	// [M-]
	// [M0]
	// [M+]
	for (int i = 0; i < this->elem.size() / 3; i++) {
		this->elem[i][0] = 1;
		this->elem[i][1] = l1[i].first;
		this->elem[i][2] = l1[i].second;
		
		this->elem[i][3] = l1[i].first * l1[i].first;
		this->elem[i][4] = l1[i].first * l1[i].second;
		this->elem[i][5] = l1[i].second * l1[i].second;

		this->elem[i][6] = l1[i].first * l1[i].first * l1[i].first;
		this->elem[i][7] = l1[i].first * l1[i].first * l1[i].second;
		this->elem[i][8] = l1[i].second * l1[i].second * l1[i].first;
		this->elem[i][9] = l1[i].second * l1[i].second * l1[i].second;
		
		this->elem[i][10] = l1[i].first * l1[i].first * l1[i].first * l1[i].first;
		this->elem[i][11] = l1[i].first * l1[i].first * l1[i].first * l1[i].second;
		this->elem[i][12] = l1[i].first * l1[i].second * l1[i].second * l1[i].first;
		this->elem[i][13] = l1[i].first * l1[i].second * l1[i].second * l1[i].second;
		this->elem[i][14] = l1[i].second * l1[i].second * l1[i].second * l1[i].second;
		
		this->elem[i][15] = l1[i].first * l1[i].first * l1[i].first * l1[i].first * l1[i].first;
		this->elem[i][16] = l1[i].first * l1[i].first * l1[i].first * l1[i].first * l1[i].second;
		this->elem[i][17] = l1[i].first * l1[i].first * l1[i].second * l1[i].second * l1[i].first;
		this->elem[i][18] = l1[i].first * l1[i].first * l1[i].second * l1[i].second * l1[i].second;
		this->elem[i][19] = l1[i].first * l1[i].second * l1[i].second * l1[i].second * l1[i].second;
		this->elem[i][20] = l1[i].second * l1[i].second * l1[i].second * l1[i].second * l1[i].second;
		
		this->elem[i][21] = l1[i].first * l1[i].first * l1[i].first * l1[i].first * l1[i].first * l1[i].first;
		this->elem[i][22] = l1[i].first * l1[i].first * l1[i].first * l1[i].first * l1[i].first * l1[i].second;
		this->elem[i][23] = l1[i].first * l1[i].first * l1[i].first * l1[i].second * l1[i].second * l1[i].first;
		this->elem[i][24] = l1[i].first * l1[i].first * l1[i].first * l1[i].second * l1[i].second * l1[i].second;
		this->elem[i][25] = l1[i].first * l1[i].first * l1[i].second * l1[i].second * l1[i].second * l1[i].second;
		this->elem[i][26] = l1[i].first * l1[i].second * l1[i].second * l1[i].second * l1[i].second * l1[i].second;
		this->elem[i][27] = l1[i].second * l1[i].second * l1[i].second * l1[i].second * l1[i].second * l1[i].second;
	}
	int j = 0;
	for (int i = this->elem.size() / 3; i < 2 * this->elem.size() / 3; i++) {
		this->elem[i][0] = 1;
		this->elem[i][1] = points[j].first;
		this->elem[i][2] = points[j].second;
		this->elem[i][3] = points[j].first * points[j].first;
		this->elem[i][4] = points[j].first * points[j].second;
		this->elem[i][5] = points[j].second * points[j].second;

		this->elem[j][6] = points[j].first * points[j].first * points[j].first;
		this->elem[j][7] = points[j].first * points[j].first * points[j].second;
		this->elem[j][8] = points[j].second * points[j].second * points[j].first;
		this->elem[j][9] = points[j].second * points[j].second * points[j].second;
		
		this->elem[i][10] = points[j].first * points[j].first * points[j].first * points[j].first;
		this->elem[i][11] = points[j].first * points[j].first * points[j].first * points[j].second;
		this->elem[i][12] = points[j].first * points[j].second * points[j].second * points[j].first;
		this->elem[i][13] = points[j].first * points[j].second * points[j].second * points[j].second;
		this->elem[i][14] = points[j].second * points[j].second * points[j].second * points[j].second;
		
		this->elem[i][15] = points[j].first * points[j].first * points[j].first * points[j].first * points[j].first;
		this->elem[i][16] = points[j].first * points[j].first * points[j].first * points[j].first * points[j].second;
		this->elem[i][17] = points[j].first * points[j].first * points[j].second * points[j].second * points[j].first;
		this->elem[i][18] = points[j].first * points[j].first * points[j].second * points[j].second * points[j].second;
		this->elem[i][19] = points[j].first * points[j].second * points[j].second * points[j].second * points[j].second;
		this->elem[i][20] = points[j].second * points[j].second * points[j].second * points[j].second * points[j].second;
		
		this->elem[i][21] = points[j].first * points[j].first * points[j].first * points[j].first * points[j].first * points[j].first;
		this->elem[i][22] = points[j].first * points[j].first * points[j].first * points[j].first * points[j].first * points[j].second;
		this->elem[i][23] = points[j].first * points[j].first * points[j].first * points[j].second * points[j].second * points[j].first;
		this->elem[i][24] = points[j].first * points[j].first * points[j].first * points[j].second * points[j].second * points[j].second;
		this->elem[i][25] = points[j].first * points[j].first * points[j].second * points[j].second * points[j].second * points[j].second;
		this->elem[i][26] = points[j].first * points[j].second * points[j].second * points[j].second * points[j].second * points[j].second;
		this->elem[i][27] = points[j].second * points[j].second * points[j].second * points[j].second * points[j].second * points[j].second;
		j++;
	}
	j = 0;
	for (int i = 2 * this->elem.size() / 3; i < this->elem.size(); i++) {
		this->elem[i][0] = 1;
		this->elem[i][1] = l3[j].first;
		this->elem[i][2] = l3[j].second;
		this->elem[i][3] = l3[j].first * l3[j].first;
		this->elem[i][4] = l3[j].first * l3[j].second;
		this->elem[i][5] = l3[j].second * l3[j].second;

		this->elem[j][6] = l3[j].first * l3[j].first * l3[j].first;
		this->elem[j][7] = l3[j].first * l3[j].first * l3[j].second;
		this->elem[j][8] = l3[j].second * l3[j].second * l3[j].first;
		this->elem[j][9] = l3[j].second * l3[j].second * l3[j].second;
		
		this->elem[i][10] = l3[j].first * l3[j].first * l3[j].first * l3[j].first;
		this->elem[i][11] = l3[j].first * l3[j].first * l3[j].first * l3[j].second;
		this->elem[i][12] = l3[j].first * l3[j].second * l3[j].second * l3[j].first;
		this->elem[i][13] = l3[j].first * l3[j].second * l3[j].second * l3[j].second;
		this->elem[i][14] = l3[j].second * l3[j].second * l3[j].second * l3[j].second;
		
		this->elem[i][15] = l3[j].first * l3[j].first * l3[j].first * l3[j].first * l3[j].first;
		this->elem[i][16] = l3[j].first * l3[j].first * l3[j].first * l3[j].first * l3[j].second;
		this->elem[i][17] = l3[j].first * l3[j].first * l3[j].second * l3[j].second * l3[j].first;
		this->elem[i][18] = l3[j].first * l3[j].first * l3[j].second * l3[j].second * l3[j].second;
		this->elem[i][19] = l3[j].first * l3[j].second * l3[j].second * l3[j].second * l3[j].second;
		this->elem[i][20] = l3[j].second * l3[j].second * l3[j].second * l3[j].second * l3[j].second;
		
		this->elem[i][21] = l3[j].first * l3[j].first * l3[j].first * l3[j].first * l3[j].first * l3[j].first;
		this->elem[i][22] = l3[j].first * l3[j].first * l3[j].first * l3[j].first * l3[j].first * l3[j].second;
		this->elem[i][23] = l3[j].first * l3[j].first * l3[j].first * l3[j].second * l3[j].second * l3[j].first;
		this->elem[i][24] = l3[j].first * l3[j].first * l3[j].first * l3[j].second * l3[j].second * l3[j].second;
		this->elem[i][25] = l3[j].first * l3[j].first * l3[j].second * l3[j].second * l3[j].second * l3[j].second;
		this->elem[i][26] = l3[j].first * l3[j].second * l3[j].second * l3[j].second * l3[j].second * l3[j].second;
		this->elem[i][27] = l3[j].second * l3[j].second * l3[j].second * l3[j].second * l3[j].second * l3[j].second;
	
		j++;
	}
}


void Matrix::fill_3L_Matrix_7th_power(const vector_pairs& points, const vector_pairs& l1, const vector_pairs& l3) {
	// заполн€ет матрицу (7 степень не€вной функции)
	// ћатрица вида 
	// [M-]
	// [M0]
	// [M+]
	for (int i = 0; i < this->elem.size() / 3; i++) {
		this->elem[i][0] = 1;
		this->elem[i][1] = l1[i].first;
		this->elem[i][2] = l1[i].second;
		
		this->elem[i][3] = l1[i].first * l1[i].first;
		this->elem[i][4] = l1[i].first * l1[i].second;
		this->elem[i][5] = l1[i].second * l1[i].second;

		this->elem[i][6] = l1[i].first * l1[i].first * l1[i].first;
		this->elem[i][7] = l1[i].first * l1[i].first * l1[i].second;
		this->elem[i][8] = l1[i].second * l1[i].second * l1[i].first;
		this->elem[i][9] = l1[i].second * l1[i].second * l1[i].second;
		
		this->elem[i][10] = l1[i].first * l1[i].first * l1[i].first * l1[i].first;
		this->elem[i][11] = l1[i].first * l1[i].first * l1[i].first * l1[i].second;
		this->elem[i][12] = l1[i].first * l1[i].second * l1[i].second * l1[i].first;
		this->elem[i][13] = l1[i].first * l1[i].second * l1[i].second * l1[i].second;
		this->elem[i][14] = l1[i].second * l1[i].second * l1[i].second * l1[i].second;
		
		this->elem[i][15] = l1[i].first * l1[i].first * l1[i].first * l1[i].first * l1[i].first;
		this->elem[i][16] = l1[i].first * l1[i].first * l1[i].first * l1[i].first * l1[i].second;
		this->elem[i][17] = l1[i].first * l1[i].first * l1[i].second * l1[i].second * l1[i].first;
		this->elem[i][18] = l1[i].first * l1[i].first * l1[i].second * l1[i].second * l1[i].second;
		this->elem[i][19] = l1[i].first * l1[i].second * l1[i].second * l1[i].second * l1[i].second;
		this->elem[i][20] = l1[i].second * l1[i].second * l1[i].second * l1[i].second * l1[i].second;
		
		this->elem[i][21] = l1[i].first * l1[i].first * l1[i].first * l1[i].first * l1[i].first * l1[i].first;
		this->elem[i][22] = l1[i].first * l1[i].first * l1[i].first * l1[i].first * l1[i].first * l1[i].second;
		this->elem[i][23] = l1[i].first * l1[i].first * l1[i].first * l1[i].second * l1[i].second * l1[i].first;
		this->elem[i][24] = l1[i].first * l1[i].first * l1[i].first * l1[i].second * l1[i].second * l1[i].second;
		this->elem[i][25] = l1[i].first * l1[i].first * l1[i].second * l1[i].second * l1[i].second * l1[i].second;
		this->elem[i][26] = l1[i].first * l1[i].second * l1[i].second * l1[i].second * l1[i].second * l1[i].second;
		this->elem[i][27] = l1[i].second * l1[i].second * l1[i].second * l1[i].second * l1[i].second * l1[i].second;
		
		this->elem[i][28] = l1[i].first * l1[i].first * l1[i].first * l1[i].first * l1[i].first * l1[i].first * l1[i].first;
		this->elem[i][29] = l1[i].first * l1[i].first * l1[i].first * l1[i].first * l1[i].first * l1[i].first * l1[i].second;
		this->elem[i][30] = l1[i].first * l1[i].first * l1[i].first * l1[i].first * l1[i].second * l1[i].second * l1[i].first;
		this->elem[i][31] = l1[i].first * l1[i].first * l1[i].first * l1[i].first * l1[i].second * l1[i].second * l1[i].second;
		this->elem[i][32] = l1[i].first * l1[i].first * l1[i].first * l1[i].second * l1[i].second * l1[i].second * l1[i].second;
		this->elem[i][33] = l1[i].first * l1[i].first * l1[i].second * l1[i].second * l1[i].second * l1[i].second * l1[i].second;
		this->elem[i][34] = l1[i].first * l1[i].second * l1[i].second * l1[i].second * l1[i].second * l1[i].second * l1[i].second;
		this->elem[i][35] = l1[i].second * l1[i].second * l1[i].second * l1[i].second * l1[i].second * l1[i].second * l1[i].second;
	}
	int j = 0;
	for (int i = this->elem.size() / 3; i < 2 * this->elem.size() / 3; i++) {
		this->elem[i][0] = 1;
		this->elem[i][1] = points[j].first;
		this->elem[i][2] = points[j].second;
		this->elem[i][3] = points[j].first * points[j].first;
		this->elem[i][4] = points[j].first * points[j].second;
		this->elem[i][5] = points[j].second * points[j].second;

		this->elem[j][6] = points[j].first * points[j].first * points[j].first;
		this->elem[j][7] = points[j].first * points[j].first * points[j].second;
		this->elem[j][8] = points[j].second * points[j].second * points[j].first;
		this->elem[j][9] = points[j].second * points[j].second * points[j].second;
		
		this->elem[i][10] = points[j].first * points[j].first * points[j].first * points[j].first;
		this->elem[i][11] = points[j].first * points[j].first * points[j].first * points[j].second;
		this->elem[i][12] = points[j].first * points[j].second * points[j].second * points[j].first;
		this->elem[i][13] = points[j].first * points[j].second * points[j].second * points[j].second;
		this->elem[i][14] = points[j].second * points[j].second * points[j].second * points[j].second;
		
		this->elem[i][15] = points[j].first * points[j].first * points[j].first * points[j].first * points[j].first;
		this->elem[i][16] = points[j].first * points[j].first * points[j].first * points[j].first * points[j].second;
		this->elem[i][17] = points[j].first * points[j].first * points[j].second * points[j].second * points[j].first;
		this->elem[i][18] = points[j].first * points[j].first * points[j].second * points[j].second * points[j].second;
		this->elem[i][19] = points[j].first * points[j].second * points[j].second * points[j].second * points[j].second;
		this->elem[i][20] = points[j].second * points[j].second * points[j].second * points[j].second * points[j].second;
		
		this->elem[i][21] = points[j].first * points[j].first * points[j].first * points[j].first * points[j].first * points[j].first;
		this->elem[i][22] = points[j].first * points[j].first * points[j].first * points[j].first * points[j].first * points[j].second;
		this->elem[i][23] = points[j].first * points[j].first * points[j].first * points[j].second * points[j].second * points[j].first;
		this->elem[i][24] = points[j].first * points[j].first * points[j].first * points[j].second * points[j].second * points[j].second;
		this->elem[i][25] = points[j].first * points[j].first * points[j].second * points[j].second * points[j].second * points[j].second;
		this->elem[i][26] = points[j].first * points[j].second * points[j].second * points[j].second * points[j].second * points[j].second;
		this->elem[i][27] = points[j].second * points[j].second * points[j].second * points[j].second * points[j].second * points[j].second;
		
		this->elem[i][28] = points[j].first * points[j].first * points[j].first * points[j].first * points[j].first * points[j].first * points[j].first;
		this->elem[i][29] = points[j].first * points[j].first * points[j].first * points[j].first * points[j].first * points[j].first * points[j].second;
		this->elem[i][30] = points[j].first * points[j].first * points[j].first * points[j].first * points[j].second * points[j].second * points[j].first;
		this->elem[i][31] = points[j].first * points[j].first * points[j].first * points[j].first * points[j].second * points[j].second * points[j].second;
		this->elem[i][32] = points[j].first * points[j].first * points[j].first * points[j].second * points[j].second * points[j].second * points[j].second;
		this->elem[i][33] = points[j].first * points[j].first * points[j].second * points[j].second * points[j].second * points[j].second * points[j].second;
		this->elem[i][34] = points[j].first * points[j].second * points[j].second * points[j].second * points[j].second * points[j].second * points[j].second;
		this->elem[i][35] = points[j].second * points[j].second * points[j].second * points[j].second * points[j].second * points[j].second * points[j].second;
		j++;
	}
	j = 0;
	for (int i = 2 * this->elem.size() / 3; i < this->elem.size(); i++) {
		this->elem[i][0] = 1;
		this->elem[i][1] = l3[j].first;
		this->elem[i][2] = l3[j].second;
		this->elem[i][3] = l3[j].first * l3[j].first;
		this->elem[i][4] = l3[j].first * l3[j].second;
		this->elem[i][5] = l3[j].second * l3[j].second;

		this->elem[j][6] = l3[j].first * l3[j].first * l3[j].first;
		this->elem[j][7] = l3[j].first * l3[j].first * l3[j].second;
		this->elem[j][8] = l3[j].second * l3[j].second * l3[j].first;
		this->elem[j][9] = l3[j].second * l3[j].second * l3[j].second;
		
		this->elem[i][10] = l3[j].first * l3[j].first * l3[j].first * l3[j].first;
		this->elem[i][11] = l3[j].first * l3[j].first * l3[j].first * l3[j].second;
		this->elem[i][12] = l3[j].first * l3[j].second * l3[j].second * l3[j].first;
		this->elem[i][13] = l3[j].first * l3[j].second * l3[j].second * l3[j].second;
		this->elem[i][14] = l3[j].second * l3[j].second * l3[j].second * l3[j].second;
		
		this->elem[i][15] = l3[j].first * l3[j].first * l3[j].first * l3[j].first * l3[j].first;
		this->elem[i][16] = l3[j].first * l3[j].first * l3[j].first * l3[j].first * l3[j].second;
		this->elem[i][17] = l3[j].first * l3[j].first * l3[j].second * l3[j].second * l3[j].first;
		this->elem[i][18] = l3[j].first * l3[j].first * l3[j].second * l3[j].second * l3[j].second;
		this->elem[i][19] = l3[j].first * l3[j].second * l3[j].second * l3[j].second * l3[j].second;
		this->elem[i][20] = l3[j].second * l3[j].second * l3[j].second * l3[j].second * l3[j].second;
		
		this->elem[i][21] = l3[j].first * l3[j].first * l3[j].first * l3[j].first * l3[j].first * l3[j].first;
		this->elem[i][22] = l3[j].first * l3[j].first * l3[j].first * l3[j].first * l3[j].first * l3[j].second;
		this->elem[i][23] = l3[j].first * l3[j].first * l3[j].first * l3[j].second * l3[j].second * l3[j].first;
		this->elem[i][24] = l3[j].first * l3[j].first * l3[j].first * l3[j].second * l3[j].second * l3[j].second;
		this->elem[i][25] = l3[j].first * l3[j].first * l3[j].second * l3[j].second * l3[j].second * l3[j].second;
		this->elem[i][26] = l3[j].first * l3[j].second * l3[j].second * l3[j].second * l3[j].second * l3[j].second;
		this->elem[i][27] = l3[j].second * l3[j].second * l3[j].second * l3[j].second * l3[j].second * l3[j].second;
		
		this->elem[i][28] = l3[j].first * l3[j].first * l3[j].first * l3[j].first * l3[j].first * l3[j].first * l3[j].first;
		this->elem[i][29] = l3[j].first * l3[j].first * l3[j].first * l3[j].first * l3[j].first * l3[j].first * l3[j].second;
		this->elem[i][30] = l3[j].first * l3[j].first * l3[j].first * l3[j].first * l3[j].second * l3[j].second * l3[j].first;
		this->elem[i][31] = l3[j].first * l3[j].first * l3[j].first * l3[j].first * l3[j].second * l3[j].second * l3[j].second;
		this->elem[i][32] = l3[j].first * l3[j].first * l3[j].first * l3[j].second * l3[j].second * l3[j].second * l3[j].second;
		this->elem[i][33] = l3[j].first * l3[j].first * l3[j].second * l3[j].second * l3[j].second * l3[j].second * l3[j].second;
		this->elem[i][34] = l3[j].first * l3[j].second * l3[j].second * l3[j].second * l3[j].second * l3[j].second * l3[j].second;
		this->elem[i][35] = l3[j].second * l3[j].second * l3[j].second * l3[j].second * l3[j].second * l3[j].second * l3[j].second;
	
		j++;
	}
}

// поиск транспонированнной матрицы
Matrix Matrix::transpose() const{
	// транспонирование матрицы
	int n = elem.size();
	int m = elem[0].size();
	Matrix N(m, n);

	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
			N[j][i] = this->elem[i][j];
	return N;
}

// поиск сопр€женной (обратной) матрицы
Matrix Matrix::inverse() {
	try {
		if (this->det() == 0) throw "Determinant = 0!";
	}
	catch (const char* exception)
	{
		cerr << "Error: " << exception << endl;
		system("pause");
		exit(EXIT_FAILURE);
	}

	int matrix_size = elem.size();
	Matrix solve_matrix(matrix_size, matrix_size * 2);

	for (int row = 0; row < matrix_size; row++) // инициализаци€ рабочей матрицы
		for (int col = 0; col < matrix_size * 2; col++)
			solve_matrix[row][col] = (col < matrix_size) ? elem[row][col] : (col == row + (matrix_size)) ? 1 : 0;

	// ѕр€мой ход
	for (int iteration = 0; iteration < matrix_size; iteration++) { // номер итерации
		if (solve_matrix[iteration][iteration] == 0) { // проверка ведущего элемента на 0
			int tmp = iteration;
			while (tmp < matrix_size) // поиск строки с ненулевым ведущим элементом
				if (solve_matrix[tmp][iteration] == 0)
					tmp++;
				else
					break;
			if (tmp == matrix_size)
				continue;
			for (int col = 0; col < 2 * matrix_size; col++) // перестановка строки
				swap(solve_matrix[tmp][col], solve_matrix[iteration][col]);
		}

		double K = solve_matrix[iteration][iteration];
		for (int col = 0; col < 2 * matrix_size; col++) // преобразование ведущей строки
			solve_matrix[iteration][col] = solve_matrix[iteration][col] / K;

		for (int row = iteration + 1; row < matrix_size; row++) { // преобразование остальных строк
			double K = solve_matrix[row][iteration] / solve_matrix[iteration][iteration];
			for (int col = 0; col < 2 * matrix_size; col++)
				solve_matrix[row][col] = solve_matrix[row][col] - (solve_matrix[iteration][col] * K);
		}
	}

	// ќбратный ход
	for (int iteration = matrix_size - 1; iteration > -1; iteration--) {

		if (solve_matrix[iteration][iteration] == 0) {
			int tmp = iteration;
			while (tmp > -1)
				if (solve_matrix[tmp][iteration] == 0)
					tmp--;
				else
					break;
			if (tmp == -1)
				continue;
			for (int col = 0; col < 2 * matrix_size; col++)
				swap(solve_matrix[tmp][col], solve_matrix[iteration][col]);
		}

		double K = solve_matrix[iteration][iteration];
		for (int col = 2 * matrix_size - 1; col > -1; col--)
			solve_matrix[iteration][col] = solve_matrix[iteration][col] / K;

		for (int row = iteration - 1; row > -1; row--)
		{
			double K = solve_matrix[row][iteration] / solve_matrix[iteration][iteration];
			for (int col = 2 * matrix_size - 1; col > -1; col--)
				solve_matrix[row][col] = solve_matrix[row][col] - (solve_matrix[iteration][col] * K);
		}
	}

	Matrix ans(matrix_size, matrix_size); // матрица с ответом
	for (int row = 0; row < matrix_size; row++) // заполнение матрицы с ответом
		for (int col = 0; col < matrix_size; col++)
			ans[row][col] = solve_matrix[row][col + matrix_size];
	return ans;
}

// многопоточное перемножение матриц
Matrix operator*(const Matrix& M, const Matrix& N) {
	int m1 = M.rows_num();
	int m2 = M.cols_num();
	int n1 = N.rows_num();
	int n2 = N.cols_num();
	if (m2 != n1) {
		cerr << "Wrong sizes!";
		exit(0);
	}
	Matrix K(m1, n2);
	int i, j, k;
	#pragma omp parallel for private(i, j, k) shared(M, N, K)
		for (i = 0; i < m1; i++) {
			for (j = 0; j < n2; j++) {
				K[i][j] = 0;
				for (k = 0; k < n1; k++) {
					K[i][j] += (M[i][k] * N[k][j]);
				}
			}
		}
	return K;
}

void Matrix::resize(int rows, int cols){
	elem.resize(rows, vector<double>(cols));
}

vector<double>& Matrix::operator[](int index){
	return elem[index];
}

const vector<double>& Matrix::operator[](int index) const{
	return elem[index];
}

Matrix mult_Matrix(const Matrix& M, const Matrix& N) {
	// перемножение двух матриц
	//M = M1;
	//N = N1;
	int m1 = M.rows_num();
	int m2 = M.cols_num();
	int n1 = N.rows_num();
	int n2 = N.cols_num();
	if (m2 != n1) {
		cerr << "Wrong sizes!";
		exit(0);
	}
	Matrix K(m1, n2);
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

Matrix& Matrix::operator=(const Matrix& M) {
	elem.resize(M.rows_num(), vector<double>(M.cols_num()));
	for (int i = 0; i < elem.size(); i++)
		for (int j = 0; j < elem[i].size(); j++)
			elem[i][j] = M[i][j];
	return *this;
}

Matrix operator-(const Matrix& M) {
	Matrix res = Matrix(M.rows_num(), M.cols_num());
	for (int i = 0; i < M.rows_num(); i++) {
		for (int j = 0; j < M.cols_num(); j++) {
			res[i][j] = -M[i][j];
		}
	}
	return res;
}

Matrix operator-(const Matrix& M1, const Matrix& M2) {
	Matrix res = Matrix(M1.rows_num(), M1.cols_num());
	for (int i = 0; i < M1.rows_num(); i++) {
		for (int j = 0; j < M1.cols_num(); j++) {
			res[i][j] = M1[i][j] - M2[i][j];
		}
	}
	return res;
}

Matrix operator+(const Matrix& M1, const Matrix& M2) {
	Matrix res = Matrix(M1.rows_num(), M1.cols_num());
	for (int i = 0; i < M1.rows_num(); i++) {
		for (int j = 0; j < M1.cols_num(); j++) {
			res[i][j] = M1[i][j] + M2[i][j];
		}
	}
	return res;
}

Matrix operator*(const double d, const Matrix& M) {
	Matrix res = Matrix(M.rows_num(), M.cols_num());
	for (int i = 0; i < M.rows_num(); i++) {
		for (int j = 0; j < M.cols_num(); j++) {
			res[i][j] = M[i][j] * 2;
		}
	}
	return res;
}

Matrix Matrix::element_power(double pow_n) {
	int n = elem.size();
	int m = elem[0].size();
	Matrix res = Matrix(n, m);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			res[i][j] = pow(elem[i][j], pow_n);
		}
	}
	return res;
}

Matrix Matrix::element_div(Matrix& m) {
	Matrix res = Matrix(rows_num(), cols_num());
	for (int i = 0; i < rows_num(); i++) {
		for (int j = 0; j < cols_num(); j++) {
			res[i][j] = elem[i][j] / m[i][j];
		}
	}
	return res;
}

Matrix Matrix::element_div(double d) {
	Matrix res = Matrix(rows_num(), cols_num());
	for (int i = 0; i < rows_num(); i++) {
		for (int j = 0; j < cols_num(); j++) {
			res[i][j] = elem[i][j] / d;;
		}
	}
	return res;
}