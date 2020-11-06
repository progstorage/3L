#include "3L.h"

class Menu;


double diff(pair<double, double> pk, pair<double, double> pk_1) {
	//производная в точке pk = (x[k], y[k]) 
	//pk_1 = (x[k+1], y[k+1])
	return ((pk_1.second - pk.second) / (pk_1.first - pk.first));
}


vector_pairs sort(const vector<double>& x, const vector<double>& y, double *d) {
	// функция сортирует исходные координаты в смысле евклидовой метрики
	// и возвращает вектор пар

	*d = 0;
	int len = x.size();

	vector_pairs sorted_points(len);
	sorted_points[0].first = x[0];
	sorted_points[0].second = y[0];
	vector<int> flags(len);

	for (int i = 0; i < len; i++) {
		flags[i] = 0;
	}
	double dist; //расстояние до первой точки
	double min_dist = sqrt((x[1] - x[0])*(x[1] - x[0]) + (y[1] - y[0])*(y[1] - y[0]));
	double max_dist = min_dist;
	int min_dist_ind = 0;
	double dydx;
	int i = 0;
	for (int k = 0; k < len - 1; k++) {
		if (flags[i] == 0) {
			for (int j = 0; j < len; j++) {
				if ((j != i) && (flags[j] == 0)) {
					dist = sqrt((x[j] - x[i])*(x[j] - x[i]) + (y[j] - y[i])*(y[j] - y[i]));
					if (dist <= min_dist) {
						min_dist = dist;
						min_dist_ind = j;
					}
					if (dist >= max_dist) { max_dist = dist; }
				}
			}
			dydx = diff(sorted_points[k], sorted_points[k + 1]);
			//cout << k << "\t(" << sorted_points[k].first << ", " << sorted_points[k].second << ")\t" << dydx << endl;
			sorted_points[k + 1].first = x[min_dist_ind];
			sorted_points[k + 1].second = y[min_dist_ind];
			flags[i] = 1;
			i = min_dist_ind;
			min_dist = 100;
		}
	}
	*d = max_dist * 0.04;
	return sorted_points;
}


void levels(vector_pairs points, double d, vector_pairs &l1, vector_pairs &l3) {
	// функция принимает на вход пары точек
	// и записывает полученные уровни в txt файл

	int len = points.size();

	double tmpx, tmpy, n;
	for (int i = 1; i < len - 1; i++) {
		tmpx = points[i + 1].second - points[i - 1].second;
		tmpy = -(points[i + 1].first - points[i - 1].first);
		n = sqrt((tmpx*tmpx + tmpy * tmpy));
		tmpx = tmpx / n;
		tmpy = tmpy / n;
		l1[i].first = points[i].first + d * tmpx;
		l1[i].second = points[i].second + d * tmpy;
		l3[i].first = points[i].first - d * tmpx;
		l3[i].second = points[i].second - d * tmpy;
		if (i == 1) {
			l1[0].first = points[0].first + d * tmpx;
			l1[0].second = points[0].second + d * tmpy;
			l3[0].first = points[0].first - d * tmpx;
			l3[0].second = points[0].second - d * tmpy;
		}
		if (i == len - 2) {
			l1[len - 1].first = points[len - 1].first + d * tmpx;
			l1[len - 1].second = points[len - 1].second + d * tmpy;
			l3[len - 1].first = points[len - 1].first - d * tmpx;
			l3[len - 1].second = points[len - 1].second - d * tmpy;
		}
	}

	ofstream outfile("out.txt", ios::out | ios::trunc);
	for (int i = 0; i < len; i++) {
		outfile << fixed << l1[i].first << "," << l1[i].second << "," << l3[i].first << "," << l3[i].second << endl;
	}
	outfile.close();
	Sleep(1000);
}


void print_Matrix(const Matrix& M) {
	// вывести матрицу
	for (const auto& v : M) {
		for (const auto& e : v)
			cout << e << "\t";
		cout << endl;
	}
	cin.get(); cin.get();
}


Matrix transpose_Matrix(const Matrix& M) {
	// транспонирование матрицы
	int n = M.size();
	int m = M[0].size();
	Matrix N(m, vector<double>(n));

	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
			N[j][i] = M[i][j];
	return N;
}


Matrix mult_Matrix(const Matrix& M, const Matrix& N) {
	// перемножение двух матриц
	int m1 = M.size();
	int m2 = M[0].size();
	int n1 = N.size();
	int n2 = N[0].size();
	if (m2 != n1) {
		cerr << "Wrong sizes!";
		exit(0);
	}
	Matrix K(m1, vector<double>(n2));

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


void fill_3L_Matrix_2nd_power(Matrix& M, const vector_pairs& points, const vector_pairs& l1, const vector_pairs& l3) {
	// заполняет матрицу (2я степень неявной функции)
	//Матрица вида 
	//[M-]
	//[M0]
	//[M+]
	for (int i = 0; i < M.size() / 3; i++) {
		M[i][0] = 1;
		M[i][1] = l1[i].first;
		M[i][2] = l1[i].second;
		M[i][3] = l1[i].first * l1[i].first;
		M[i][4] = l1[i].first * l1[i].second;
		M[i][5] = l1[i].second * l1[i].second;
	}
	int j = 0;
	for (int i = M.size() / 3; i < 2 * M.size() / 3; i++) {
		M[i][0] = 1;
		M[i][1] = points[j].first;
		M[i][2] = points[j].second;
		M[i][3] = points[j].first * points[j].first;
		M[i][4] = points[j].first * points[j].second;
		M[i][5] = points[j].second * points[j].second;
		j++;
	}
	j = 0;
	for (int i = 2 * M.size() / 3; i < M.size(); i++) {
		M[i][0] = 1;
		M[i][1] = l3[j].first;
		M[i][2] = l3[j].second;
		M[i][3] = l3[j].first * l3[j].first;
		M[i][4] = l3[j].first * l3[j].second;
		M[i][5] = l3[j].second * l3[j].second;
		j++;
	}
}


double det_Matrix(const Matrix& M) {
	// определитель Метод Гаусса
	const double EPS = 1E-9;
	int n = M.size();
	Matrix N(M);	// Копирование матрицы
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


Matrix inverse_Matrix(const Matrix& M) {
	// находит обратную матрицу
	try
	{
		if (det_Matrix(M) == 0) throw "Determinant of a matrix = 0";
	}
	catch (const char* exception)
	{
		cerr << "Error: " << exception << endl;
		system("pause");
		exit(EXIT_FAILURE);
	}

	int matrix_size = M.size();
	Matrix solve_matrix(matrix_size, vector<double>(matrix_size * 2));

	for (int row = 0; row < matrix_size; row++) // инициализация рабочей матрицы
		for (int col = 0; col < matrix_size * 2; col++)
			solve_matrix[row][col] = (col < matrix_size) ? M[row][col] : (col == row + (matrix_size)) ? 1 : 0;

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

	Matrix ans(matrix_size, vector<double>(matrix_size)); // матрица с ответом
	for (int row = 0; row < matrix_size; row++) // заполнение матрицы с ответом
		for (int col = 0; col < matrix_size; col++)
			ans[row][col] = solve_matrix[row][col + matrix_size];
	return ans;
}


void solve_system(Matrix& M, double d) {
	// Решает систему уравнений - находит коэффициенты многочлена и записывает их в файл "coef.txt"
	double len = M.size();
	Matrix b(1, vector<double>(len));
	for (int i = 0; i < len / 3; i++) { b[0][i] = d; }
	for (int i = len / 3; i < 2 * len / 3; i++) { b[0][i] = 0; }
	for (int i = 2 * len / 3; i < len; i++) { b[0][i] = -d; }
	Matrix a = mult_Matrix(mult_Matrix(inverse_Matrix(mult_Matrix(transpose_Matrix(M), M)), transpose_Matrix(M)), transpose_Matrix(b));

	ofstream outfile("coef.txt", ios::out | ios::trunc);
	for (int i = 0; i < 6; i++) {
		outfile << fixed << a[i][0] << endl;
	}
	outfile.close();
	Sleep(1000);
}


Menu::Menu(string path_to_dir) {
	dir = path_to_dir;
	get_files(path_to_dir);
	start();
	to_draw = true;
}


void Menu::get_files(string path_to_dir) {
	for (const auto & f : fs::directory_iterator(path_to_dir)) {
		files.push_back(f.path().string());
	}
}


void Menu::print_menu() {
	cout << "Directory: " << dir << endl;
	cout << "File" << "\t\t\t" << "N" << endl;
	
	int i;
	for (i = 0; i < files.size(); i++) {
		cout << files[i] << "\t" << i + 1 << endl;
	}

	cout << "Exit" << "\t\t\t" << i + 1 << endl;
}


void Menu::select_file() {
	cout << endl << "Enter file number:\t";
	cin >> n;
	if (n <= files.size()) {
		test_file = files[n - 1];
	}
	else if (n == files.size() + 1) {
		to_draw = false;
	}
	else {
		cout << endl << "Error! Try again!";
		select_file();
	}
	cout << "\n-------------------------\n\n";
}


void Menu::start() {
	print_menu();
	select_file();
	read_file();
}


void Menu::split(string line, double* x, double* y) {
	auto pos = line.find(";");
	if (pos != string::npos) {
		*x = stod(line.substr(0, pos));
		*y = stod(line.substr(pos + 1));
	}
}


void Menu::read_file() {
	string line;
	double tmpx;
	double tmpy;
	ifstream in(test_file);		//Чтение данных из файла test_file
	if (in.is_open()) {
		while (getline(in, line)) {
			split(line, &tmpx, &tmpy);
			x.push_back(tmpx);
			y.push_back(tmpy);
		}
	}
	in.close();
}
