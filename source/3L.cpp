#include "3L.h"

class Menu;


double diff(pair<double, double> pk, pair<double, double> pk_1) {
	//производная в точке pk = (x[k], y[k]) 
	//pk_1 = (x[k+1], y[k+1])
	return ((pk_1.second - pk.second) / (pk_1.first - pk.first));
}


int old_sort(const vector<double>& x, const vector<double>& y, double* d) {
	// функция сортирует исходные координаты в смысле евклидовой метрики
	// и возвращает вектор пар

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
			sorted_points[k + 1].first = x[min_dist_ind];
			sorted_points[k + 1].second = y[min_dist_ind];
			flags[i] = 1;
			i = min_dist_ind;
			min_dist = 100;
		}
	}
	*d = max_dist * 0.04;

	points_clusters_array.push_back(sorted_points);//отсортированная группа точек добавляется в массив групп
	points_clusters_count++;

	ofstream outfile("sorted_points/sorted_points.txt", ios::out | ios::trunc);
	for (int i = 0; i < len; i++) {
		outfile << fixed << sorted_points[i].first << "," << sorted_points[i].second << endl;
	}

	return 0;
}


int sort_method_2(const vector<double>& x, const vector<double>& y, double* d) {
	// функция сортирует исходные координаты в смысле евклидовой метрики
	// если расстояние до следующей точки "слишком большое" (в данном случае - больше расстояния до первой точки),
	// то все точки разбиваются на отдельные группы, в которых расстояния от каждой точки до следующей не превышает этого значения
	// => формируется глобальный вектор отсортированных векторов пар, а также каждая группа точек сохраняется в файл {число}.txt

	*d = 0;
	int len = x.size();

	vector_pairs sorted_points(len);
	sorted_points[0].first = x[0];
	sorted_points[0].second = y[0];
	vector<int> flags(len);

	for (int i = 0; i < len; i++) {
		flags[i] = 0;
	}
	double dist, dist_to_fitrs;
	double min_dist = sqrt((x[1] - x[0])*(x[1] - x[0]) + (y[1] - y[0])*(y[1] - y[0]));
	double max_dist = min_dist;
	int min_dist_ind = 0;
	int i = 0;
	int retflag;
	string type = ".txt";
	string dir = "sorted_points/";
	for (int k = 0; k < len - 1; k++) {
		if (flags[i] == 0) {
			flags[i] = 1;
			for (int j = 0; j < len; j++) {
				if ((j != i) && (flags[j] == 0)) {
					retflag = 1;//если нашлась точка, в которой еще не были
					dist = sqrt((x[j] - x[i])*(x[j] - x[i]) + (y[j] - y[i])*(y[j] - y[i]));
					if (dist <= min_dist) {
						min_dist = dist;
						min_dist_ind = j;
					}
					if (dist >= max_dist) { max_dist = dist; }
				}
			}
			dist_to_fitrs = sqrt((x[min_dist_ind] - x[0])*(x[min_dist_ind] - x[0]) + (y[min_dist_ind] - y[0])*(y[min_dist_ind] - y[0]));
			if (min_dist <= dist_to_fitrs) {
				sorted_points[k + 1].first = x[min_dist_ind]; 
				sorted_points[k + 1].second = y[min_dist_ind];
				i = min_dist_ind;
				min_dist = 100;
				retflag = 0;
			}
			else {//если расстояние до следующей точки больше, чем расстояние до первой точки (вообще говоря, если больше диаметра группы,
				  //но как определить диаметр группы пока не понятно
				vector_pairs tmp_sort_points(k + 1);
				for (int j = 0; j < k + 1; j++) {
					tmp_sort_points[j] = sorted_points[j];
				}
				points_clusters_array.push_back(tmp_sort_points);//отсортированная группа точек добавляется в массив групп
				points_clusters_count++;
				auto name = dir + to_string(points_clusters_count) + type;
				ofstream outfile(name, ios::out | ios::trunc);
				for (int i = 0; i < k + 1; i++) {
					outfile << fixed << sorted_points[i].first << "," << sorted_points[i].second << endl;
				}

				flags[i] = 1;
				flags[min_dist_ind] = 1;
				vector<double> new_x;//оставшиеся точки формируют новые массивы координат х и у
				vector<double> new_y;
				for (int i = 0; i < len; i++) {
					if (flags[i] == 0) {
						new_x.push_back(x[i]);
						new_y.push_back(y[i]);
					}
				}

				if (new_x.size() != 0) { sort_method_2(new_x, new_y, d); }//рекурсивный вызов функции для оставшихся точек
				else { return 0; }
			}
		}
	}
	if (retflag != 1) {//если перебраны все точки
		points_clusters_array.push_back(sorted_points);
		points_clusters_count++;
		auto name = dir + to_string(points_clusters_count) + type;
		ofstream outfile(name, ios::out | ios::trunc);
		for (int i = 0; i < sorted_points.size(); i++) {
			outfile << fixed << sorted_points[i].first << "," << sorted_points[i].second << endl;
		}
		return 0;
	}
	*d = max_dist * 0.04;
	return 0;
}



void levels(vector_pairs points, double d, vector_pairs& l1, vector_pairs& l3, int num) {
	// функция принимает на вход пары точек
	// и записывает полученные уровни в txt файл
	string type = ".txt";
	string dir = "levels/";

	int len = points.size();

	double tmpx, tmpy, n;
	for (int i = 1; i < len - 1; i++) {
		tmpx = points[i + 1].second - points[i - 1].second;
		tmpy = -(points[i + 1].first - points[i - 1].first);
		n = sqrt((tmpx*tmpx + tmpy*tmpy));
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

	auto name = dir + to_string(num) + type;
	ofstream outfile(name, ios::out | ios::trunc);
	for (int i = 0; i < len; i++) {
		outfile << fixed << l1[i].first << "," << l1[i].second << "," << l3[i].first << "," << l3[i].second << endl;
	}
	outfile.close();
}


void print_Matrix(const Matrix& M) {
	// вывести матрицу
	for (const auto& v : M) {
		for (const auto& e : v)
			cout << e << "\t";
		cout << endl;
	}
	//cin.get(); cin.get();
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
	//M = M1;
	//N = N1;
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


Matrix generate_random_Matrix(const int n, const int m) {
	// генерирует случайную матрицу размера n x m
	Matrix res(n, vector<double>(m));
	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
			res[i][j] = (rand() % 10000 - 5000) / 100;
	return res;
}


Matrix mult_Matrix_multithread(const Matrix& M, const Matrix& N) {
	int m1 = M.size();
	int m2 = M[0].size();
	int n1 = N.size();
	int n2 = N[0].size();
	if (m2 != n1) {
		cerr << "Wrong sizes!";
		exit(0);
	}
	Matrix K(m1, vector<double>(n2));
	
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


void solve_system(Matrix& M, double d, int num) {
	// Решает систему уравнений - находит коэффициенты многочлена и записывает их в файл "coef.txt"
	string type = ".txt";
	string dir = "coefs/";
	double len = M.size();
	Matrix b(1, vector<double>(len));
	for (int i = 0; i < len / 3; i++) { b[0][i] = d; }
	for (int i = len / 3; i < 2 * len / 3; i++) { b[0][i] = 0; }
	for (int i = 2 * len / 3; i < len; i++) { b[0][i] = -d; }

	Matrix a;

	#if TEST_MODE == 1
		clock_t start_1, end_1, start_2, end_2;
	
		// многопоточное умножение матрицы на вектор
		start_1 = clock();
		a = mult_Matrix_multithread(mult_Matrix_multithread(inverse_Matrix(mult_Matrix_multithread(transpose_Matrix(M), M)), transpose_Matrix(M)), transpose_Matrix(b));		
		end_1 = clock() - start_1;
			
		// однопоточное
		start_2 = clock();
		Matrix _a = mult_Matrix(mult_Matrix(inverse_Matrix(mult_Matrix(transpose_Matrix(M), M)), transpose_Matrix(M)), transpose_Matrix(b));
		end_2 = clock() - start_2;
		
		cout << "\nSingle-threaded multiplication execution time:\t" << end_2 / (double)CLOCKS_PER_SEC
			 << "\nMulti-threaded multiplication execution time:\t" << end_1 / (double)CLOCKS_PER_SEC;
	
		if (_a == a)
			cout << "\nBoth methods gave the same result!" << endl;
		else
			cout << "\nDifferent results, something went wrong!" << endl;
	#elif TEST_MODE == 2
		// big-size matrix testing
		clock_t start_1, end_1, start_2, end_2;

		// многопоточное умножение матрицы на вектор
		start_1 = clock();
		a = mult_Matrix_multithread(generate_random_Matrix(size_1, size_2), generate_random_Matrix(size_2, size_3));
		end_1 = clock() - start_1;

		// однопоточное
		start_2 = clock();
		Matrix _a = mult_Matrix(generate_random_Matrix(size_1, size_2), generate_random_Matrix(size_2, size_3));
		end_2 = clock() - start_2;

		cout << "\nSingle-threaded multiplication execution time:\t" << end_2 / (double)CLOCKS_PER_SEC
			 << "\nMulti-threaded multiplication execution time:\t" << end_1 / (double)CLOCKS_PER_SEC;
		cout << endl;
	#else
		a = mult_Matrix_multithread(mult_Matrix_multithread(inverse_Matrix(mult_Matrix_multithread(transpose_Matrix(M), M)), transpose_Matrix(M)), transpose_Matrix(b));
	#endif

	auto name = dir + to_string(num) + type;
	ofstream outfile(name, ios::out | ios::trunc);
	for (int i = 0; i < 6; i++) {
		outfile << fixed << a[i][0] << endl;
	}
	outfile.close();
}


void fill_levels_multithread(double d) {
	int i;
	vector_pairs l1, l3;
	Matrix M;

	#if TEST_MODE
		clock_t start_1, end_1, start_2, end_2;
		start_1 = clock();

		// многопоточное создание уровней
		#pragma omp parallel for private(i) shared(points_clusters_array, l1, l3)
			for (i = 0; i < points_clusters_count; i++) {
				l1.resize(points_clusters_array[i].size());
				l3.resize(points_clusters_array[i].size());

				levels(points_clusters_array[i], d, l1, l3, i + 1);

				M.resize(points_clusters_array[i].size() * 3, vector<double>(6));
				fill_3L_Matrix_2nd_power(M, points_clusters_array[i], l1, l3);
				solve_system(M, d, i + 1);

				l1.clear();
				l3.clear();
			}
		end_1 = clock() - start_1;

		// однопоточное
		start_2 = clock();
		for (int j = 0; j < points_clusters_count; j++) {
			l1.resize(points_clusters_array[j].size());
			l3.resize(points_clusters_array[j].size());
			levels(points_clusters_array[j], d, l1, l3, j + 1);

			M.resize(points_clusters_array[j].size() * 3, vector<double>(6));
			fill_3L_Matrix_2nd_power(M, points_clusters_array[j], l1, l3);
			solve_system(M, d, j + 1);

			l1.clear();
			l3.clear();
		}
		end_2 = clock() - start_2;

		cout << "\nSingle-threaded level creation time:\t" << end_2 / (double)CLOCKS_PER_SEC
			 << "\nMulti-threaded level creation time:\t" << end_1 / (double)CLOCKS_PER_SEC;

	#else
		#pragma omp parallel for private(i) shared(points_clusters_array, l1, l3)
			for (i = 0; i < points_clusters_count; i++) {
				l1.resize(points_clusters_array[i].size());
				l3.resize(points_clusters_array[i].size());

				levels(points_clusters_array[i], d, l1, l3, i + 1);

				M.resize(points_clusters_array[i].size() * 3, vector<double>(6));
				fill_3L_Matrix_2nd_power(M, points_clusters_array[i], l1, l3);
				solve_system(M, d, i + 1);

				l1.clear();
				l3.clear();
			}
	#endif
}

void fill_3L_Matrix_2nd_power(Matrix& M, const vector_pairs& points, const vector_pairs& l1, const vector_pairs& l3) {
	// заполняет матрицу (2я степень неявной функции)
	// Матрица вида 
	// [M-]
	// [M0]
	// [M+]
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
		exit(0);
	}
	else {
		cout << endl << "Error! Try again!";
		select_file();
	}
	cout << "\n-------------------------\n\n";
}

void Menu::select_execution_method() {
	cout << endl << "Enter 1 - old method, 2 - new method:\t";
	cin >> execution_method;
}

void Menu::start() {
	print_menu();
	select_file();
	read_file();
	select_execution_method();
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
