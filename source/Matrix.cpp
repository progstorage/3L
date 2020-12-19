#include "Matrix.h"

// конструктор
Matrix::Matrix(int rows, int cols) {
	// this->elem = new std::vector<std::vector<double> >(rows, std::vector<double> (cols, 0));
	this->elem.resize(rows, vector<double>(cols));
}

Matrix::Matrix(Matrix& M) {
	int rows = M.rows_num();
	int cols = M.cols_num();
	// this->elem = new std::vector<std::vector<double> >(rows, std::vector<double> (cols));	
	this->elem.resize(rows, vector<double>(cols));
	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
			elem[i][j] = M[i][j];
}

// деструктор
Matrix::~Matrix() {
	// for (int i = 0; i < elem.size(); i++)
	// 	delete elem[i];
	// delete elem;
}

// возврат размерности
int Matrix::rows_num() const { return elem.size(); }
int Matrix::cols_num() const { return elem[0].size(); }

// поиск определителя
double Matrix::det() {
	// определитель Метод Гаусса
	const double EPS = 1E-9;
	int n = elem.size();
	Matrix N(*this);	// Копирование матрицы
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
		//swap(N[i], N[k]);
		if (i != k)
			det = -det;
		det *= N[i][i];
		for (int j = i + 1; j < n; ++j) {
			N[i][j] /= N[i][i];
		}
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

// заполнение матрицы
void Matrix::fill_3L_Matrix_2nd_power(const vector_pairs& points, const vector_pairs& l1, const vector_pairs& l3) {
	// заполняет матрицу (2я степень неявной функции)
	// Матрица вида 
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

// поиск транспонированнной матрицы
Matrix Matrix::transpose() {
	// транспонирование матрицы
	int n = elem.size();
	int m = elem[0].size();
	Matrix N(m, n);

	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
			N[j][i] = this->elem[i][j];
	return N;
}

// поиск сопряженной (обратной) матрицы
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

	for (int row = 0; row < matrix_size; row++) // инициализация рабочей матрицы
		for (int col = 0; col < matrix_size * 2; col++)
			solve_matrix[row][col] = (col < matrix_size) ? elem[row][col] : (col == row + (matrix_size)) ? 1 : 0;

	// Прямой ход
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

	// Обратный ход
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
Matrix& Matrix::operator*(Matrix& N) const {
	int m1 = elem.size();
	int m2 = elem[0].size();
	int n1 = N.rows_num();
	int n2 = N.cols_num();
	if (m2 != n1) {
		cerr << "Wrong sizes!";
		exit(0);
	}
	Matrix K(m1, n2);
	int i, j, k;
#pragma omp parallel for private(i, j, k) shared(elem, N, K)
	for (i = 0; i < m1; i++) {
		for (j = 0; j < n2; j++) {
			K[i][j] = 0;
			for (k = 0; k < n1; k++) {
				K[i][j] += (elem[i][k] * N[k][j]);
			}
		}
	}
	return K;
}

void Matrix::resize(int rows, int cols) {
	elem.resize(rows, vector<double>(cols));
}

vector<double>& Matrix::operator[](int index) {
	return elem[index];
}

Matrix mult_Matrix(Matrix& M, Matrix& N) {
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